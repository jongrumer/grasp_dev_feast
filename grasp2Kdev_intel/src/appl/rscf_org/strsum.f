************************************************************************
      SUBROUTINE STRSUM
      IMPLICIT REAL*8          (A-H,O-Z)

*   Generates the first part of  rscf92.sum  (on stream 24).
*
*   Call(s) to: [LIB92] CALEN, CONVRT.
*
*   Written by Farid A. Parpia            Last revision: 26 Sep 1993
*
************************************************************************

      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)

      POINTER (PNTRWT,RWTDUMMY)
      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL DIAG,LFORDR
      CHARACTER*256 RECORD
      CHARACTER*26 CDATA, ctime*8, cdate*8
      CHARACTER*2 CLEVEL,NH
*
      POINTER (PWEIGH,WEIGHT(1))
      POINTER (PCCMIN,ICCMIN(1))
      POINTER (PNTRPF,PF(NNNP,1))
      POINTER (PNTRQF,QF(NNNP,1))
*
      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /DEF2/C
     :      /DEF5/PNTRWT,PWEIGH
     :      /DEF7/PCCMIN,NCMIN,NCMAX
     :      /FOPARM/ICCUT
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /MCPB/DIAG,LFORDR
     :      /NPAR/PARM(2),NPARM
     :      /NPOT/ZZ(NNNP),NNUC
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
     :      /WFAC/WFACT
*
*   Both the nuclear charge and the number of electrons are
*   known at this point; load IONCTY with the ionicity
*
      IONCTY = NINT (Z) - NELEC
*
*   Get the date and time of day; make this information the
*   header of the summary file
*
!      CALL CALEN (CTIME, CDATE)
!      WRITE (24,*) 'RSCF92 run at ',CTIME,' on ',CDATE,'.'
*
*   Write out the basic dimensions of the electron cloud
*
      WRITE (24,*)
      CALL CONVRT (NELEC, RECORD, LENTH)
      WRITE (24,*) 'There are '//RECORD(1:LENTH)
     :           //' electrons in the cloud'
      CALL CONVRT (NCF, RECORD, LENTH)
      WRITE (24,*) ' in '//RECORD(1:LENTH)
     :           //' relativistic CSFs'
      CALL CONVRT (NW, RECORD, LENTH)
      WRITE (24,*) ' based on '//RECORD(1:LENTH)
     :           //' relativistic subshells.'
*
*   If the CSFs are not treated uniformly, write out an
*   informative message
*
      IF (LFORDR) THEN
         WRITE (24,*)
         CALL CONVRT (ICCUT, RECORD, LENTH)
         WRITE (24,*) ' CSFs 1--'//RECORD(1:LENTH)//' constitute'
     :              //' the zero-order space;'
      ENDIF
*
*   Write out the nuclear parameters
*
      WRITE (24,*)
      WRITE (24,300) Z
      IF (EMN .EQ. 0.D0) THEN
         WRITE (24,*) ' the nucleus is stationary;'
      ELSE
         WRITE (24,301) EMN
      ENDIF
      IF (NPARM .EQ. 2) THEN
         WRITE (24,*) ' Fermi nucleus:'
         WRITE (24,302) PARM(1), PARM(2)
         CALL CONVRT (NNUC, RECORD, LENTH)
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
      IF (HP .EQ. 0.D0) THEN
         WRITE (24,305) RNT, H, N
      ELSE
         WRITE (24,306) RNT, H, HP, N
      ENDIF
      WRITE (24,307) R(1), R(2), R(N)
      WRITE (24,*)
*
*  (E)AL calculation, returns here
*
      IF (NCMIN .EQ. 0) THEN
         WRITE (24,*) '(E)AL calculation.'
         RETURN
      ENDIF
*
*  Info exclusively for EOL calculations
*
      IF (NCMIN .EQ. 1) THEN
         WRITE (24,*) 'OL calculation.'
         CALL CONVRT (ICCMIN(1), RECORD, LENTH)
         WRITE (24,*) 'Level '//RECORD(1:LENTH)//' will be optimised.'
      ELSE
         WRITE (24,*) 'EOL calculation.'
         CALL CONVRT (NCMIN, RECORD, LENTH)
         WRITE (24,*) RECORD(1:LENTH)//' levels will be optimised;'
         RECORD (1:20) = ' their indices are: '
         IEND = 20
         DO 2 I = 1, NCMIN
            IBEG = IEND + 1
            CALL CONVRT (ICCMIN(I), CLEVEL, LENTH)
            IF (I .NE. NCMIN) THEN
               IEND = IBEG + LENTH + 1
               RECORD(IBEG:IEND) = CLEVEL(1:LENTH)//', '
            ELSE
               IEND = IBEG+LENTH
               RECORD(IBEG:IEND) = CLEVEL(1:LENTH)//'.'
            ENDIF
            IF (IEND .GE. 120) THEN
               WRITE (24,*) RECORD(1:IEND)
               RECORD(1:2) = '  '
               IEND = 2
            ENDIF
    2    CONTINUE
         IF (IEND .NE. 2) WRITE (24,*) RECORD(1:IEND)
         IF (WEIGHT(1) .EQ. -1.D0) THEN
            WRITE (24,*) 'Each is assigned its statistical weight;'
         ELSEIF (WEIGHT(1) .EQ. -2.D0) THEN
            WRITE (24,*) 'All levels are weighted equally;'
         ELSE
            WRITE (24,*) ' weighted as follows:'
            WRITE (24,*) (WEIGHT(I), I = 1, NCMIN)
         ENDIF
      ENDIF

  300 FORMAT ('The atomic number is ',1F14.10,';')
  301 FORMAT (' the mass of the nucleus is ',1PD19.12,
     :        ' electron masses;')
  302 FORMAT ('  c =',1P,1D19.12,' Bohr radii,'
     :       /'  a =',   1D19.12,' Bohr radii;')
  303 FORMAT ('Speed of light = ',3PD19.12,' atomic units.')
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

      RETURN
      END
