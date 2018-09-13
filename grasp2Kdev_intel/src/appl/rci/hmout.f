************************************************************************
*                                                                      *
      SUBROUTINE HMOUT (IMCDF)
*                                                                      *
*   Routine for printing the Hamiltonian matrix.  File IMCDF must be   *
*   positioned correctly before a call is made to this module.         *
*                                                                      *
*   Call(s) to: [LIB92]: CONVRT.                                       *
*                                                                      *
*   Written by Farid A Parpia             Last revision: 14 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)

      integer*8 nelmnt

Cww      INTEGER PIENDC,PNTRIQ
      POINTER (PIENDC,ENDCDUMMY)
      POINTER (PNTRIQ,RIQDUMMY)
      CHARACTER*8 CIR,CIC
*
      POINTER (PNTEMT,EMT(*))
      POINTER (PNIROW,IROW(*))
*
      COMMON/HMAT/PNTEMT,PIENDC,PNIROW,NELMNT
     :      /ORB2/NCF,NW,PNTRIQ
*
      DO 2 IC = 1,NCF
         READ(IMCDF) NELC,ELSTO,(EMT(IR),IR = 1,NELC),
     :                          (IROW(IR),IR = 1,NELC)
         DO 1 IR = 1,NELC
            CALL CONVRT (IROW(IR),CIR,LENTHR)
            CALL CONVRT (     IC ,CIC,LENTHC)
            WRITE (99,300) CIR(1:LENTHR),CIC(1:LENTHC),EMT(IR)
    1    CONTINUE
    2 CONTINUE
*
  300 FORMAT (' H (',A,',',A,') = ',1P,1D19.12)
*
      END
