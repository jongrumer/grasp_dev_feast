************************************************************************
*                                                                      *
      SUBROUTINE CONNECT
*                                                                      *
*   The position of an orbital in the merged list is connected to      *
*   the positions in the initial and final state lists                 *
*                                                                      *
*   Written by Per Jonsson                                             *
*                                                                      *
************************************************************************

      IMPLICIT DOUBLEPRECISION(A-H,O-Z)
*
      include 'parameters.def'
CFF      PARAMETER (NNNW = 120)

Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)

      COMMON/ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB2/NCF,NW,PNTRIQ
      COMMON/ORB4II/NPII(NNNW),NAKII(NNNW)
     :      /ORB2II/NCFII,NWII
     :      /ORBCONII/NNII(NNNW)
      COMMON/ORB4FF/NPFF(NNNW),NAKFF(NNNW)
     :      /ORB2FF/NCFFF,NWFF
     :      /ORBCONFF/NNFF(NNNW)
*
*     Initialize
*
      DO I = 1,NW
        NNII(I) = 0.D0
        NNFF(I) = 0.D0
      ENDDO
*
*   Loop over the orbitals in the merged list
*
      DO I=1,NW
        DO J=1,NWII
          IF (NP(I).EQ.NPII(J).AND.NAK(I).EQ.NAKII(J)) NNII(I)=J
        ENDDO
        DO J=1,NWFF
          IF (NP(I).EQ.NPFF(J).AND.NAK(I).EQ.NAKFF(J)) NNFF(I)=J
        ENDDO
      ENDDO

      RETURN
      END
