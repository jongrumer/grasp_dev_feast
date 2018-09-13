!***********************************************************************
!                                                                      *
      INTEGER FUNCTION ITJPO (ICSF) 
!                                                                      *
!   ITJPO is the value of 2J+1 for CSF number ICSF.                    *
!                                                                      *
!   Call(s) to: [LIB92]: IUNPCK.                                       *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 02 Nov 1992   *
!   Modified by G. Gaigalas                                 May 2011   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:48:45   2/14/04  
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE STAT_C,      ONLY: JCUPA
      USE IOUNIT_C,    ONLY: ISTDE
      USE orb_C,       ONLY: NCF
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER       :: ICSF 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: NNNW = 120 
!-----------------------------------------------
      IF (ICSF>=1 .AND. ICSF<=NCF) THEN 
!        bit extraction will not give the proper sign.
!GG         ITJPO = IBITS(JCUPA((NNNW - 1)/4 + 1,ICSF),8*MOD(NNNW - 1,4),8) 
        itjpo = jcupa(NNNW,icsf)
        IF (ITJPO > 127) ITJPO = 256 - ITJPO 
        ITJPO = IABS (ITJPO)
!        ITJPO = IUNPCK (JCUPA(1,ICSF),NNNW)
!        IF (itjpob .ne. itjpo) then
!          print *, jcupa((nnnw-1)/4+1,ICSF)
!          print *, 'ITJPO', nnnw,icsf,itjpo,itjpob
!          stop
!        END IF
      ELSE 
         WRITE (ISTDE, *) 'ITJPO: Argument ICSF is out of range.' 
         STOP  
      ENDIF 
!
      RETURN  
      END FUNCTION ITJPO 
