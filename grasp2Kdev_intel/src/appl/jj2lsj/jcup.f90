!***********************************************************************
!                                                                      *
      INTEGER FUNCTION JCUP (LOC, ICSF) 
!                                                                      *
!   JCUP is the 2J+1 value of the LOCth nontrivial intermediate ang-   *
!   ular momentum in CSF  ICSF.                                        *
!                                                                      *
!   Call(s) to: [LIB92]: IUNPCK.                                       *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 02 Nov 1992   *
!   Modified by G. Gaigalas                                 May 2011   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  08:16:51   2/21/04  
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE IOUNIT_C,   ONLY: ISTDE
      use orb_C,      ONLY: NW, NCF
      use stat_C,     ONLY: JCUPA
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: LOC 
      INTEGER              :: ICSF 
!-----------------------------------------------
!
! Packed array declared I*4 and a common (iounit) added
! XHH 1997.02.12
!
      IF (LOC>=1 .AND. LOC<=NW-1) THEN 
         IF (ICSF>=1 .AND. ICSF<=NCF) THEN 
!cff        JCUP = IUNPCK (JCUPA(1,ICSF),LOC)
!GG            JCUP = IBITS(JCUPA((LOC - 1)/4 + 1,ICSF),8*MOD(LOC - 1,4),8) 
            jcup = jcupa(loc,icsf)
         ELSE 
            WRITE (ISTDE, *) 'JCUP: Argument ICSF is out of range.' 
            STOP  
         ENDIF 
      ELSE 
         WRITE (ISTDE, *) 'JCUP: Argument LOC is out of range.' 
         STOP  
      ENDIF 
!
      RETURN  
      END FUNCTION JCUP 
