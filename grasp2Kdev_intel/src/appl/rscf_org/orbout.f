************************************************************************

      SUBROUTINE ORBOUT (rwffile2)
      IMPLICIT REAL*8          (A-H, O-Z)
      CHARACTER rwffile2*(*)

*   This routine is modified to:
*   Open, write and close the RWF file. This routine is called after
*   the end of each iteration. So the intermediate (iteration) result 
*   will most probably be safe if the execution is interrupted.
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 18 Dec 1992   *
*   Modified by Xinghong He               Last revision: 05 Aug 1998   *
*                                                                      *
************************************************************************

      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)

      POINTER (PNTRIQ,RIQDUM)
      POINTER (PNTRPF,PF(NNNP,1))
      POINTER (PNTRQF,QF(NNNP,1))

      COMMON/GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /ORB1/E(NNNW),GAMA(NNNW)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
!-----------------------------------------------------------------------
      OPEN (23, FILE=rwffile2, STATUS='UNKNOWN', FORM='UNFORMATTED')
      WRITE (23) 'G92RWF'
      DO J = 1, NW
         MFJ = MF(J)
         WRITE (23) NP(J), NAK(J), E(J), MFJ
         WRITE (23) PZ(J), (PF(I,J), I=1,MFJ), (QF(I,J), I=1,MFJ)
         WRITE (23) (R(I), I=1,MFJ)   ! This is a waste of resources
      ENDDO
      CLOSE (23)

      RETURN
      END
