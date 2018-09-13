************************************************************************
*                                                                      *
      SUBROUTINE STRSUM
*                                                                      *
*   Generates the first part of  erwf.sum  (on stream 24).             *
*                                                                      *
*   Call(s) to: [LIB92] CALEN, CONVRT.                                 *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 17 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      CHARACTER*256 RECORD
      CHARACTER*26 CDATA
*
      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /DEF2/C
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /NPAR/PARM(2),NPARM
     :      /NPOT/ZZ(NNNP),NNUC
     :      /ORB2/NCF,NW,PNTRIQ
*
*   Get the date and time of day; make this information the
*   header of the summary file
*
      CALL CALEN (CTIME,CDATE)
      WRITE (24,*) 'ERWF run at ',CTIME,' on ',CDATE,'.'
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
*   Write out the nuclear parameters
*
      WRITE (24,*)
      WRITE (24,300) Z
      IF (NPARM .EQ. 2) THEN
         WRITE (24,*) ' Fermi nucleus:'
         WRITE (24,302) PARM(1),PARM(2)
         CALL CONVRT (NNUC,RECORD,LENTH)
         WRITE (24,*) ' there are '//RECORD(1:LENTH)
     :              //' tabulation points in the nucleus.'
      ELSE
         WRITE (24,*) ' point nucleus.'
      ENDIF
*
*   Write out the physical effects specifications
*
      WRITE (24,*)
      WRITE (24,303) C
*
*   Write out the parameters of the radial grid
*
      WRITE (24,*)
      IF (HP .EQ. 0.0D 00) THEN
         WRITE (24,305) RNT,H,N
      ELSE
         WRITE (24,306) RNT,H,HP,N
      ENDIF
      WRITE (24,307) R(1),R(2),R(N)
*
      WRITE (24,*)
*
      RETURN
*
  300 FORMAT ('The atomic number is ',1F14.10,';')
  302 FORMAT ('  c =',1P,1D19.12,' Bohr radii,'
     :       /'  a =',   1D19.12,' Bohr radii;')
  303 FORMAT ('Speed of light = ',1PD19.12,' atomic units.')
  305 FORMAT ( 'Radial grid: R(I) = RNT*(exp((I-1)*H)-1),',
     :         ' I = 1, ..., N;'
     :       //' RNT  = ',1P,D19.12,' Bohr radii;'
     :        /' H    = ',   D19.12,' Bohr radii;'
     :        /' N    = ',1I4,';')
  306 FORMAT ( 'Radial grid: ln(R(I)/RNT+1)+(H/HP)*R(I) = (I-1)*H,',
     :         ' I = 1, ..., N;'
     :       //' RNT  = ',1P,D19.12,' Bohr radii;'
     :        /' H    = ',   D19.12,' Bohr radii;'
     :        /' HP   = ',   D19.12,' Bohr radii;'
     :        /' N    = ',1I4,';')
  307 FORMAT ( ' R(1) = ',1P,1D19.12,' Bohr radii;'
     :        /' R(2) = ',   1D19.12,' Bohr radii;'
     :        /' R(N) = ',   1D19.12,' Bohr radii.')
*
      END
