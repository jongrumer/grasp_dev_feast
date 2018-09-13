************************************************************************
*                                                                      *
      FUNCTION ITJPO (ICSF)
*                                                                      *
*   ITJPO is the value of 2J+1 for CSF number ICSF.                    *
*                                                                      *
*   Call(s) to: [LIB92]: IUNPCK.                                       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 02 Nov 1992   *
*                                                                      *
************************************************************************
*
! Packed array JCUPA declared I*4 and common/iounit/ added
! XHH 1997.02.12
Cww      INTEGER PNTJQS,PNTRIQ
      POINTER (PNTJQS,JQSDUMMY)
      POINTER (PNTRIQ,RIQDUMMY)
*
*
      include 'parameters.def'
CGG      INTEGER NNNW
CGG      PARAMETER (NNNW = 120)
CGG      INTEGER NNNWP
CGG      PARAMETER (NNNWP = 30)

      INTEGER*4 JCUPA
      POINTER (PNJCUP,JCUPA(NNNWP,*))

      COMMON/ORB2/NCF,NW,PNTRIQ
     :      /STAT/PNTJQS,PNJCUP
      COMMON/iounit/istdi,istdo,istde
*
      IF ((ICSF .GE. 1) .AND. (ICSF .LE. NCF)) THEN
*        bit extraction will not give the proper sign.
         ITJPO = IBITS(JCUPA((NNNW-1)/4+1,ICSF),8*MOD(NNNW-1,4),8)
         IF (ITJPO .gt. 127) ITJPO = 256-ITJPO
*        ITJPO = IUNPCK (JCUPA(1,ICSF),NNNW)
*        ITJPO = ABS (ITJPO)
*        IF (itjpob .ne. itjpo) then
*          print *, jcupa((nnnw-1)/4+1,ICSF)
*          print *, 'ITJPO', nnnw,icsf,itjpo,itjpob
*          stop
*        END IF
      ELSE
         WRITE(istde,*) 'ITJPO: Argument ICSF is out of range.'
         STOP
      ENDIF
*
      RETURN
      END
