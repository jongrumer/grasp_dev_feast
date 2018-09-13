************************************************************************
*                                                                      *
      SUBROUTINE RKINTC (IA,IB,IC,ID,K,TEGRAL)
*                                                                      *
*                         k                                            *
*   This routine returns R (abcd) integrals.                           *
*                                                                      *
*   Written by Per Jonsson                                             *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      PARAMETER (KMAX = 20)
      LOGICAL FOUND,FIRST
Cww        INTEGER PNTRIQ
      POINTER(PNTRIQ,RIQDUMMY(*))                 
*
      POINTER (PCTEILRK,INDTEIRK(*))
      POINTER (PCTEVLRK,VALTEIRK(*))
*
      COMMON/CTEILSRK/PCTEILRK,PCTEVLRK
     :      /ORB2/NCF,NW,PNTRIQ
     :      /KKSTART/KSTART(0:KMAX)
*
C      PARAMETER (KEY = 121)
      KEY = NW + 1
*
*   Ensure that the indices are in `canonical' order
*   Compute the composite (packed) index
*


      IF (IA.GT.IC) THEN
        ISWAP = IC
        IC = IA
        IA = ISWAP
      ENDIF
      IF (IB.GT.ID) THEN
        ISWAP = ID
        ID = IB
        IB = ISWAP
      ENDIF
      IF (IA.GT.IB) THEN
        ISWAP = IB
        IB = IA
        IA = ISWAP
        ISWAP = ID
        ID = IC
        IC = ISWAP
      ENDIF

      INDEX = ((IA*KEY+IB)*KEY+IC)*KEY+ID
* 
      JL = KSTART(K)
      JU = KSTART(K+1) - 1

      IF (INDEX.LT.INDTEIRK(JL).OR.INDEX.GT.INDTEIRK(JU)) THEN
        WRITE(*,*) 'Something wrong in rkintc'
        STOP
      ENDIF
*
*   The index is within the range of the indices stored; search
*   for it in the list of indices
*
    1 IF (JU-JL .GT. 1) THEN
        JM = (JU+JL)/2
        IF (INDTEIRK(JM) .GT. INDEX) THEN
          JU = JM
        ELSE
          JL = JM
        ENDIF
        GOTO 1
      ENDIF
*
*   The range is bracketed to the extent possible
*
      IF (INDEX .EQ. INDTEIRK(JU)) THEN
        LOC = JU
      ELSEIF (INDEX .EQ. INDTEIRK(JL)) THEN
        LOC = JL
      ELSE
        WRITE(*,*) K,IA,IB,IC,ID,INDEX
        WRITE(*,*) 'Rkintc Integral not found'
        STOP
      ENDIF
*
*   Return the value of the integral
*   from storage

      TEGRAL = VALTEIRK(LOC)
*
      RETURN
      END
