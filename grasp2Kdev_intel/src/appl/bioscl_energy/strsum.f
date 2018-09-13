************************************************************************
*                                                                      *
      SUBROUTINE STRSUM(NAME,INPCI,ILBL)
*                                                                      *
*   Generates the first part of  oscl92.sum  (on stream 24).           *
*                                                                      *
*   Call(s) to: [LIB92]: CALEN, CONVRT, WGHTD5.                        *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 28 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)

      include 'parameters.def'
CFF      PARAMETER (NNNP = 590) 
CFF      PARAMETER (NNN1 = NNNP+10)
CFF      PARAMETER (NNNW = 120)

      LOGICAL LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
      CHARACTER*24 NAME(2)
      CHARACTER*15 RECORD 
      CHARACTER*8 CTIME,CDATE
      CHARACTER*2 NHII,NHFF
*
*  Initial state pointers
*
      POINTER (PNEVALII,EVALII(1))
      POINTER (PNEVECII,EVECII(1))
      POINTER (PNIVECII,IVECII(1))
      POINTER (PNTRPFII,PFII(NNNP,1)),(PNTRQFII,QFII(NNNP,1))
      POINTER (PIATJPII,IATJPOII(1)),(PIASPAII,IASPARII(1))
*
*  Final state pointers
*
      POINTER (PNEVALFF,EVALFF(1))
      POINTER (PNEVECFF,EVECFF(1))
      POINTER (PNIVECFF,IVECFF(1))
      POINTER (PNTRPFFF,PFFF(NNNP,1)),(PNTRQFFF,QFFF(NNNP,1))
      POINTER (PIATJPFF,IATJPOFF(1)),(PIASPAFF,IASPARFF(1))
*
      COMMON/DECIDE/LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
     :      /DEF1/EMN,IONCTY,NELEC,Z
     :      /DEF2/C
     :      /FOPARM/ICCUT
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /NPAR/PARM(2),NPARM
*
*  Common relevant for the initial state
*
      COMMON/EIGVALII/EAVII,PNEVALII
     :      /EIGVECII/PNEVECII
     :      /DEF1II/EMNII,IONCTYII,NELECII,ZII
     :      /ORB1II/EII(NNNW),GAMAII(NNNW)
     :      /ORB2II/NCFII,NWII
     :      /ORB4II/NPII(NNNW),NAKII(NNNW)
     :      /ORB10II/NHII(NNNW)
     :      /PRNTII/NVECII,PNIVECII,NVECMXII
     :      /SYMAII/PIATJPII,PIASPAII
     :      /WAVEII/PZII(NNNW),PNTRPFII,PNTRQFII,MFII(NNNW)
*
*  Common relevant for the final state
*
      COMMON/EIGVALFF/EAVFF,PNEVALFF
     :      /EIGVECFF/PNEVECFF
     :      /DEF1FF/EMNFF,IONCTYFF,NELECFF,ZFF
     :      /ORB1FF/EFF(NNNW),GAMAFF(NNNW)
     :      /ORB2FF/NCFFF,NWFF
     :      /ORB4FF/NPFF(NNNW),NAKFF(NNNW)
     :      /ORB10FF/NHFF(NNNW)
     :      /PRNTFF/NVECFF,PNIVECFF,NVECMXFF
     :      /SYMAFF/PIATJPFF,PIASPAFF
     :      /WAVEFF/PZFF(NNNW),PNTRPFFF,PNTRQFFF,MFFF(NNNW)

CGG   lbl beginning
      REAL Lev_ENER_1,Lev_ENER_2
      CHARACTER*1 Lev_par_1,Lev_par_2
      CHARACTER*4 Lev_J_1,Lev_J_2
      CHARACTER*64 string_CSF1,string_CSF2
      POINTER(PNELev_POS_1,Lev_POS_1(1))
      POINTER(PNELev_J_1,Lev_J_1(1))
      POINTER(PNELev_Par_1,Lev_Par_1(1))
      POINTER(PNELev_ENER_1,Lev_ENER_1(1))
      POINTER(PNEstring_CSF1,string_CSF1(1))
      COMMON/JJ2LSJ1/ NVECTOTI,IOPEN_STATUS1,PNELev_POS_1,PNELev_J_1,
     :                PNELev_Par_1,PNELev_ENER_1,PNEstring_CSF1
*
      POINTER(PNELev_POS_2,Lev_POS_2(1))
      POINTER(PNELev_J_2,Lev_J_2(1))
      POINTER(PNELev_Par_2,Lev_Par_2(1))
      POINTER(PNELev_ENER_2,Lev_ENER_2(1))
      POINTER(PNEstring_CSF2,string_CSF2(1))
      COMMON/JJ2LSJ2/ NVECTOTF,IOPEN_STATUS2,PNELev_POS_2,PNELev_J_2,
     :                PNELev_Par_2,PNELev_ENER_2,PNEstring_CSF2
CGG   lbl end
      I = INDEX(NAME(1),' ')
      J = INDEX(NAME(2),' ')
      IF(ILBL . EQ. 0) THEN
         IF (INPCI.EQ.0) THEN
           OPEN(UNIT=24,FILE=NAME(1)(1:I-1)//'.'//NAME(2)(1:J-1)//'.ct',
     :          FORM='FORMATTED',STATUS='UNKNOWN')
         ELSE
           OPEN(UNIT=24,FILE=NAME(1)(1:I-1)//'.'//NAME(2)(1:J-1)//'.t',
     :          FORM='FORMATTED',STATUS='UNKNOWN')
         ENDIF
      ELSE IF(ILBL . EQ. 1) THEN
         IF(IOPEN_STATUS1.EQ.0 .AND. IOPEN_STATUS2 .EQ.0) THEN
            IF (INPCI.EQ.0) THEN
               OPEN(UNIT=32,
     :              FILE=NAME(1)(1:I-1)//'.'//NAME(2)(1:J-1)//'.ct.lsj',
     :              FORM='FORMATTED',STATUS='UNKNOWN')
            ELSE
               OPEN(UNIT=32,
     :              FILE=NAME(1)(1:I-1)//'.'//NAME(2)(1:J-1)//'.t.lsj',
     :              FORM='FORMATTED',STATUS='UNKNOWN')
            ENDIF
         END IF
      END IF
*
*   Get the date and time of day; make this information the
*   header of the summary file
*
Cww      CALL CALEN (CTIME,CDATE)
Cww      WRITE (24,*) 'OSCL92 run at ',CTIME,' on ',CDATE,'.'
*
*   Write out the basic dimensions of the initial state electron cloud
*
Cww      WRITE (24,*)
Cww      CALL CONVRT (NELECII,RECORD,LENTH)
Cww      WRITE (24,*) 'There are '//RECORD(1:LENTH)
Cww     :           //' electrons in the initial state cloud'
Cww      CALL CONVRT (NCFII,RECORD,LENTH)
Cww      WRITE (24,*) ' in '//RECORD(1:LENTH)
Cww     :           //' relativistic CSFs'
Cww      CALL CONVRT (NWII,RECORD,LENTH)
Cww      WRITE (24,*) ' based on '//RECORD(1:LENTH)
Cww     :           //' relativistic subshells.'
*
*   Write out the basic dimensions of the final state electron cloud
*
Cww      WRITE (24,*)
Cww      CALL CONVRT (NELECFF,RECORD,LENTH)
Cww      WRITE (24,*) 'There are '//RECORD(1:LENTH)
Cww     :           //' electrons in the final state cloud'
Cww      CALL CONVRT (NCFFF,RECORD,LENTH)
Cww      WRITE (24,*) ' in '//RECORD(1:LENTH)
Cww     :           //' relativistic CSFs'
Cww      CALL CONVRT (NWFF,RECORD,LENTH)
Cww      WRITE (24,*) ' based on '//RECORD(1:LENTH)
Cww     :           //' relativistic subshells.'
*
*   If the CSFs are not treated uniformly, write out an
*   informative message
*
Cww      IF (LFORDR) THEN
Cww         WRITE (24,*)
Cww         CALL CONVRT (ICCUT,RECORD,LENTH)
Cww         WRITE (24,*) ' CSFs 1--'//RECORD(1:LENTH)//' constitute'
Cww     :              //' the zero-order space;'
Cww      ENDIF
*
*   Write out the nuclear parameters
*
Cww      WRITE (24,*)
Cww      WRITE (24,300) Z
Cww      IF (NPARM .EQ. 2) THEN
Cww         WRITE (24,*) 'Fermi nucleus:'
Cww         WRITE (24,301) PARM(1),PARM(2)
Cww      ELSE
Cww         WRITE (24,*) ' point nucleus.'
Cww      ENDIF
*
*   Write out the physical effects specifications
*
Cww      WRITE (24,*)
Cww      WRITE (24,305) C
*
*   Write out the parameters of the radial grid
*
Cww      WRITE (24,*)
Cww      IF (HP .EQ. 0.0D 00) THEN
Cww         WRITE (24,306) RNT,H,N
Cww      ELSE
Cww         WRITE (24,307) RNT,H,HP,N
Cww      ENDIF
Cww      WRITE (24,308) R(1),R(2),R(N)
*
*   Write out the orbital properties
*
Cww      WRITE (24,*)
Cww      WRITE (24,*) 'Initial state subshell radial wavefunction summary:'
Cww      WRITE (24,*)
Cww      WRITE (24,309)
Cww      WRITE (24,*)
Cww      DO 1 I = 1,NWII
Cww         WRITE (24,310) NPII(I),NHII(I),EII(I),PZII(I),
Cww     :                  GAMAII(I),PFII(2,I),QFII(2,I),MFII(I)
Cww    1 CONTINUE
*
Cww      WRITE (24,*)
Cww      WRITE (24,*) 'Final state subshell radial wavefunction summary:'
Cww      WRITE (24,*)
Cww      WRITE (24,309)
Cww      WRITE (24,*)
Cww      DO 2 I = 1,NWFF
Cww         WRITE (24,310) NPFF(I),NHFF(I),EFF(I),PZFF(I),
Cww     :                  GAMAFF(I),PFFF(2,I),QFFF(2,I),MFFF(I)
Cww    2 CONTINUE
*
*   Write the list of eigenpair indices for the initial state
*
c     WRITE (24,*)
c     CALL ENGOUT1 (EAVII,EVALII,IATJPOII,IASPARII,IVECII,NVECII,3,1)
*
*   Write the list of eigenpair indices for the final state
*
c     WRITE (24,*)
c     CALL ENGOUT1 (EAVFF,EVALFF,IATJPOFF,IASPARFF,IVECFF,NVECFF,3,2)
*
      RETURN
*
  300 FORMAT ('The atomic number is ',1F14.10,';')
  301 FORMAT ('  c =',1P,1D19.12,' Bohr radii,'
     :       /'  a =',   1D19.12,' Bohr radii;')
  305 FORMAT ('Speed of light = ',1PD19.12,' atomic units.')
  306 FORMAT ( 'Radial grid: R(I) = RNT*(exp((I-1)*H)-1),',
     :         ' I = 1, ..., N;'
     :       //' RNT  = ',1P,D19.12,' Bohr radii;'
     :        /' H    = ',   D19.12,' Bohr radii;'
     :        /' N    = ',1I4,';')
  307 FORMAT ( 'Radial grid: ln(R(I)/RNT+1)+(H/HP)*R(I) = (I-1)*H,',
     :         ' I = 1, ..., N;'
     :       //' RNT  = ',1P,D19.12,' Bohr radii;'
     :        /' H    = ',   D19.12,' Bohr radii;'
     :        /' HP   = ',   D19.12,' Bohr radii;'
     :        /' N    = ',1I4,';')
  308 FORMAT ( ' R(1) = ',1P,1D19.12,' Bohr radii;'
     :        /' R(2) = ',   1D19.12,' Bohr radii;'
     :        /' R(N) = ',   1D19.12,' Bohr radii.')
  309 FORMAT (' Subshell',11X,'e',20X,'p0',18X,
     :        'gamma',19X,'P(2)',18X,'Q(2)',10X,'MTP')
  310 FORMAT (3X,1I2,1A2,1X,1P,5(3X,1D19.12),3X,1I3)
*
      END
