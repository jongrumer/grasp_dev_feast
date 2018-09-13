************************************************************************
*                                                                      *
      SUBROUTINE TCSL(N)
*                                                                      *
*   This subroutine transfers data to the initial and final state      *
*   common blocks                                                      * 
*                                                                      *
*   Written by Per Jonsson                Last revision: June   1996   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)

      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)

      CHARACTER*2 NH,NHII,NHFF

      POINTER (PNTRIQ,PNTDUMMY)

      COMMON/ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
     :      /ORB10/NH(NNNW)
*
*  Initial state commons
*
      COMMON/ORB2II/NCFII,NWII
     :      /ORB4II/NPII(NNNW),NAKII(NNNW)
     :      /ORB5II/NKLII(NNNW),NKJII(NNNW)
     :      /ORB10II/NHII(NNNW)
*
*  Final state commons
*
      COMMON/ORB2FF/NCFFF,NWFF
     :      /ORB4FF/NPFF(NNNW),NAKFF(NNNW)
     :      /ORB5FF/NKLFF(NNNW),NKJFF(NNNW)
     :      /ORB10FF/NHFF(NNNW)
*
      IF (N.EQ.1) THEN
         NCFII = NCF
         NWII = NW
         DO I = 1,NW
            NPII(I) = NP(I)
            NAKII(I) = NAK(I)
            NKLII(I) = NKL(I)
            NKJII(I) = NKJ(I)
            NHII(I) = NH(I)
         ENDDO
      ELSE
         NCFFF = NCF
         NWFF = NW
         DO I = 1,NW
            NPFF(I) = NP(I)
            NAKFF(I) = NAK(I)
            NKLFF(I) = NKL(I)
            NKJFF(I) = NKJ(I)
            NHFF(I) = NH(I)
         ENDDO
      ENDIF

      RETURN
      END 
