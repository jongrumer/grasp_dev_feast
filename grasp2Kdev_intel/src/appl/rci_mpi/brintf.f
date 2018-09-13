************************************************************************
*                                                                      *
      FUNCTION BRINTF (ITYPE,IA,IB,IC,ID,K)
*                                                                      *
*   Computes integrals for the transverse photon interaction.          *
*                                                                      *
*   Call(s) to: [RCI92]: BESSEL, BRRA.                                 *
*                                                                      *
*                                           Last update: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
*
      COMMON/STOR/KEEP(2,2)
*
      include 'parameters.def'
CGG      PARAMETER (KEYORB = 121)
      PARAMETER (KEY = KEYORB, KEY2 = KEY*KEY)
*
      GOTO (1,4,3,6,2,5),ITYPE
*
*   Type 1 and 5 integrals require j(k), n(k) Bessel fuctions
*   Type 5 integrals only require w = wab Bessel functions
*
    1 IF ( ((IA .NE. IC) .OR. (IB .NE. ID)) .OR.
     :     ((IA .NE. ID) .OR. (IC .EQ. IB))) CALL BESSEL (IC,ID,1,2,K)
    2 CALL BESSEL (IA,IB,1,1,K)
      GOTO 6
*
*   Type 3 integrals require j(k), n(k) Bessel functions for either
*   w = wab or w = cd whichever is non-zero.
*
    3 IF (IA .NE. IB) CALL BESSEL (IA,IB,1,1,K)
      IF (IC .NE. ID) CALL BESSEL (IC,ID,1,2,K)
      GOTO 6
*
*   Type 2 and 6 integrals require j(k), n(k) and j(k+2), n(k+2)
*   Bessel fuctions
*   Type 6 integrals only require w = wab Bessel functions.
*
    4 CONTINUE
      IF (((IA .NE. IC) .OR. (IB .NE. ID)) .OR.
     :    ((IA .NE. ID) .OR. (IC .NE. IB))) THEN
*
         ICOD = MAX (IC,ID)+KEY*MIN (IC,ID)
         ICOD1 = ICOD+KEY2*(K-1)
         ICOD2 = ICOD+KEY2*(K+1)
         ICOD = MAX (IA,IB)+KEY*MIN (IA,IB)
         ICOD3 = ICOD+KEY2*(K-1)
         ICOD4 = ICOD+KEY2*(K+1)
         IF ((ICOD1 .EQ. KEEP(1,2)) .AND.
     :       (ICOD2 .EQ. KEEP(2,2)) .AND.
     :       (ICOD3 .EQ. KEEP(1,1)) .AND.
     :       (ICOD4 .EQ. KEEP(2,1))) GOTO 6
         IF ((ICOD1 .EQ. KEEP(1,1)) .AND.
     :       (ICOD2 .EQ. KEEP(2,1)) .AND.
     :       (ICOD3 .EQ. KEEP(1,2)) .AND.
     :       (ICOD4 .EQ. KEEP(2,2))) GOTO 6
         CALL BESSEL (IC,ID,1,2,K-1)
         CALL BESSEL (IC,ID,2,2,K+1)
      ENDIF
*
    5 CALL BESSEL (IA,IB,1,1,K-1)
      CALL BESSEL (IA,IB,2,1,K+1)
*
*   Compute the integral
*
    6 BRINTF = BRRA (ITYPE,IA,IB,IC,ID,K)
*
      RETURN
      END
