************************************************************************
*                                                                      *
      SUBROUTINE STRSUM
*                                                                      *
*   Generates the first part of  sms92.sum  (on stream 24).            *
*                                                                      *
*   Call(s) to: [LIB92] CALEN, CONVRT, WGHTD5.                         *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 28 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = 600)
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
      CHARACTER*256 RECORD
      CHARACTER*26 CDATA
      CHARACTER*2 NH
*
      POINTER (PNEVAL,EVAL(*))
      POINTER (PNEVEC,EVEC(*))
      POINTER (PNIVEC,IVEC(*))
      POINTER (PNTRPF,PF(NNNP,*)),(PNTRQF,QF(NNNP,*))
      POINTER (PIATJP,IATJPO(*)),(PIASPA,IASPAR(*))
*
      COMMON/DECIDE/LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
     :      /DEF1/EMN,IONCTY,NELEC,Z
     :      /DEF2/C
     :      /EIGVAL/EAV,PNEVAL
     :      /EIGVEC/PNEVEC
     :      /FOPARM/ICCUT
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /NPAR/PARM(2),NPARM
     :      /NPOT/ZZ(NNNP),NNUC
     :      /NSMDAT/SMSI,SMSD,SMSQ
     :      /ORB1/E(NNNW),GAMA(NNNW)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /PRNT/NVEC,PNIVEC,NVECMX
     :      /SYMA/PIATJP,PIASPA
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
*
*   Get the date and time of day; make this information the
*   header of the summary file
*
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
      IF (LFORDR) THEN
         WRITE (24,*)
         CALL CONVRT (ICCUT,RECORD,LENTH)
         WRITE (24,*) ' CSFs 1--'//RECORD(1:LENTH)//' constitute'
     :              //' the zero-order space;'
      ENDIF
*
*   Write out the nuclear parameters
*
      WRITE (24,*)
      WRITE (24,300) Z
      IF (NPARM .EQ. 2) THEN
         WRITE (24,*) 'Fermi nucleus:'
         WRITE (24,301) PARM(1),PARM(2)
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
      WRITE (24,305) C
*
*   Write out the parameters of the radial grid
*
      WRITE (24,*)
      IF (HP .EQ. 0.0D 00) THEN
         WRITE (24,306) RNT,H,N
      ELSE
         WRITE (24,307) RNT,H,HP,N
      ENDIF
      WRITE (24,308) R(1),R(2),R(N)
*
*   Write out the orbital properties
*
      WRITE (24,*)
      WRITE (24,*) 'Subshell radial wavefunction summary:'
      WRITE (24,*)
      WRITE (24,309)
      WRITE (24,*)
      DO 1 I = 1,NW
         WRITE (24,310) NP(I),NH(I),E(I),PZ(I),
     :                  GAMA(I),PF(2,I),QF(2,I),MF(I)
    1 CONTINUE
*
*   Write the list of eigenpair indices
*
      WRITE (24,*)
      CALL ENGOUT (EAV,EVAL,IATJPO,IASPAR,IVEC,NVEC,3)
      CALL WGHTD5
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
