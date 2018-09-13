!***********************************************************************
!                                                                      *
      INTEGER FUNCTION ISPAR (ICSF) 
!                                                                      *
!   ISPAR is the value of P for CSF number ICSF.                       *
!                                                                      *
!   Call(s) to: [LIB92]: IUNPCK.                                       *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 02 Nov 1992   *
!   Modified by G. Gaigalas                                 May 2011   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:48:41   2/14/04  
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE STAT_C,      ONLY: JCUPA
      USE IOUNIT_C,    ONLY: ISTDE 
      USE ORB_C,       ONLY: NCF
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
!
! The packed integer array jcupa declared as I*4
! COMMON/iounit/ added
! XHH 1997.02.12
!
      IF (ICSF>=1 .AND. ICSF<=NCF) THEN 
!cff     ISPAR = IUNPCK (JCUPA(1,ICSF),NNNW)
!        .. note this bit extraction does not preserve sign
!GG         ISPAR = IBITS(JCUPA((NNNW - 1)/4 + 1,ICSF),8*MOD(NNNW - 1,4),8) 
         ispar = jcupa(NNNW,icsf)
         IF (ISPAR > 127) ISPAR = ISPAR - 256 
         ISPAR = SIGN(1,ISPAR) 
      ELSE 
         WRITE (ISTDE, *) 'ISPAR: Argument ICSF is out of range.' 
         STOP  
      ENDIF 
!
      RETURN  
      END FUNCTION ISPAR 
