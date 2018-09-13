************************************************************************
*                                                                      *
      SUBROUTINE STRSUM
*                                                                      *
*   Generates the first part of the  .sum  file on stream 24.          *
*                                                                      *
*   Call(s) to: [LIB92] CALEN, CONVRT.                                 *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 19 Nov 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H,O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL DIAG,LFORDR
      CHARACTER*256 RECORD
      CHARACTER*26 CDATA
      CHARACTER*2 NH
*
      POINTER (pncfblk, ncfblk(0:*))
      COMMON/hblock/nblock, pncfblk

      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /FOPARM/ICCUT(100)
     :      /MCPB/DIAG,LFORDR
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
*
*   Get the date and time of day; make this information the
*   header of the summary file
*
      CALL CALEN (CTIME,CDATE)                                          
      WRITE (24,*) 'GENMCP run at ',CTIME,' on ',CDATE,'.'
*
*   Write out the basic dimensions of the electron cloud
*
      WRITE (24,*)
      CALL CONVRT (NELEC,RECORD,LENTH)
      WRITE (24,*) 'There are '//RECORD(1:LENTH)
     :           //' electrons in the cloud'
      CALL CONVRT (NCF,RECORD,LENTH)
      WRITE (24,*) ' in '//RECORD(1:LENTH)
     :           //' relativistic CSFs'
      CALL CONVRT (NW,RECORD,LENTH)
      WRITE (24,*) ' based on '//RECORD(1:LENTH)
     :           //' relativistic subshells.'
*
*   If the CSFs are not treated uniformly, write out an
*   informative message
*
      IF (DIAG) THEN
         WRITE (24,*)
         WRITE (24,*) 'Only diagonal matrix elements are computed.'
      ELSE
         IF (LFORDR) THEN
            do i = 1,nblock
               WRITE (24,*)
               CALL CONVRT (ICCUT(i),RECORD,LENTH)
               WRITE (24,*) 'CSFs 1--'//RECORD(1:LENTH)//' constitute'
     :                 //' the zero-order space.'
            end do
         ENDIF
      ENDIF
*
      RETURN
      END
