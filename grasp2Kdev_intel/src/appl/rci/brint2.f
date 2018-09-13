************************************************************************
*                                                                      *
      SUBROUTINE BRINT2 (IA,IB,IC,ID,K,TEGRAL)
*                                                                      *
*   Returns integrals for the transverse photon interaction.           *
*                                                                      *
*   Written by Per Jonsson                   Octaober 2014             *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      PARAMETER (KMAX = 20)
Cww   INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)

      POINTER (PINDT1,INDTP1(*))
      POINTER (PVALT1,VALTP1(*))
      POINTER (PINDT2,INDTP2(*))
      POINTER (PVALT2,VALTP2(*))
      POINTER (PINDT3,INDTP3(*))
      POINTER (PVALT3,VALTP3(*))
      POINTER (PINDT4,INDTP4(*))
      POINTER (PVALT4,VALTP4(*))
      POINTER (PINDT5,INDTP5(*))
      POINTER (PVALT5,VALTP5(*))
      POINTER (PINDT6,INDTP6(*))
      POINTER (PVALT6,VALTP6(*))
*
      COMMON/BILST/PINDT1,PINDT2,PINDT3,PINDT4,PINDT5,PINDT6,
     :             PVALT1,PVALT2,PVALT3,PVALT4,PVALT5,PVALT6,
     :             NDTPA(6),NTPI(6),FIRST(6)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /BALLOC/NB(2)
      COMMON/KKSTARTBREIT/KSTARTBREIT1(0:KMAX),KSTARTBREIT2(0:KMAX)
*
      KEY = NW + 1

      INDEX = ((IA*KEY+IB)*KEY+IC)*KEY+ID
*
      JL = KSTARTBREIT2(K)
      JU = KSTARTBREIT2(K+1) - 1

      IF (INDEX.LT.INDTP2(JL).OR.INDEX.GT.INDTP2(JU)) THEN
        WRITE(*,*) 'Something wrong in brint2'
        STOP
      ENDIF
*
*   The index is within the range of the indices stored; search
*   for it in the list of indices
*
    1 IF (JU-JL .GT. 1) THEN
        JM = (JU+JL)/2
        IF (INDTP2(JM) .GT. INDEX) THEN
          JU = JM
        ELSE
          JL = JM
        ENDIF
        GOTO 1
      ENDIF
*
*   The range is bracketed to the extent possible
*
      IF (INDEX .EQ. INDTP2(JU)) THEN
        LOC = JU
      ELSEIF (INDEX .EQ. INDTP2(JL)) THEN
        LOC = JL
      ELSE
        WRITE(*,*) K,IA,IB,IC,ID,INDEX
        WRITE(*,*) 'Brint2 Integral not found'
        STOP
      ENDIF
*
*   Return the value of the integral
*   from storage

      TEGRAL = VALTP2(LOC)
*
      RETURN
      END

