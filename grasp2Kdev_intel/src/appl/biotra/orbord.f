************************************************************************
*                                                                      *
      SUBROUTINE ORBORD(N)
*                                                                      *
*   This subroutine checks the orbital order. Normal and reversed      *
*   order can be used. If reveresed order JA and JB must be            *
*   permuted in the subroutine ti1tv                                   *
*                                                                      *
*   Written by Per Jonsson                Last revision: Feb    1997   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)

      include 'parameters.def
CGG      PARAMETER (NNNW = 120)

      CHARACTER*2 NH

      POINTER (PNTRIQ,PNTDUMMY)

      COMMON/ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
     :      /ORB10/NH(NNNW)
     :      /ORB5FF/NKLFF(NNNW),NKJFF(NNNW)
     :      /ORB10FF/NHFF(NNNW)
*

      RETURN
      END 
