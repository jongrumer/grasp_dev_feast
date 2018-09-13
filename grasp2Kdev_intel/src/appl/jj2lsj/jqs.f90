!!***********************************************************************
!                                                                      *
      INTEGER FUNCTION JQS (IWHICH, ISUBSH, ICSF) 
!                                                                      *
!   JQS is a subshell quantum number for subshell ISUBSH in configu-   *
!   ration state function  ICSF:  the seniority if IWHICH is 1;  the   *
!   quantum number w if IWHICH is 2, and 2J+1 if IWHICH is 3.          *
!                                                                      *
!   Call(s) to: [LIB92]: IUNPCK.                                       *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 02 Nov 1992   *
!   Modified by G. Gaigalas                                 May 2011   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:48:50   2/14/04  
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE STAT_C,  ONLY: JQSA
!GG      USE IOUNIT_C 
!GG      use orb_C, ONLY: NCF, NW
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER             :: IWHICH 
      INTEGER, INTENT(IN) :: ISUBSH 
      INTEGER             :: ICSF 
!-----------------------------------------------
!
! Packed array declared I*4 and a common added (iounit)
! XHH 1997.02.12
!
!
!     IF ((IWHICH .GE. 1) .AND. (IWHICH .LE. 3)) THEN
!        IF ((ISUBSH .GE. 1) .AND. (ISUBSH .LE. NW)) THEN
!           IF ((ICSF .GE. 1) .AND. (ICSF .LE. NCF)) THEN
!cff           JQS = IUNPCK (JQSA(1,IWHICH,ICSF),ISUBSH)
!GG      JQS = IBITS(JQSA((ISUBSH - 1)/4 + 1,IWHICH,ICSF),8*MOD(ISUBSH - 1,4),8) 
      jqs = jqsa(isubsh,iwhich,icsf)
!           ELSE
!              WRITE(istde,*) 'JQS: Argument ICSF is out of range.'
!              STOP
!           ENDIF
!        ELSE
!           WRITE(istde,*) 'JQS: Argument ISUBSH is out of range.'
!           STOP
!        ENDIF
!     ELSE
!        WRITE(istde,*) 'JQS: Argument IWHICH is out of range.'
!        STOP
!     ENDIF
!
      RETURN  
      END FUNCTION JQS 
