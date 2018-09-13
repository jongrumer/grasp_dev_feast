************************************************************************
*                                                                      *
      SUBROUTINE ENDSUM
*                                                                      *
*   Generates the last part of  rscf92.sum  (on stream 24).            *
*                                                                      *
*   Call(s) to: [LIB92]: ENGOUT, RINT.                                 *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 24 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H,O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      CHARACTER*2 NH
*
      POINTER (PCCMIN,ICCMIN(1))
      POINTER (PNEVAL,EVAL(1))
      POINTER (PIATJP,IATJPO(1))
      POINTER (PIASPA,IASPAR(1))
      POINTER (PNTRPF,PF(NNNP,1))
      POINTER (PNTRQF,QF(NNNP,1))
*
      COMMON/DEF7/PCCMIN,NCMIN,NCMAX
     :      /EIGVAL/EAV,PNEVAL
     :      /ORB1/E(NNNW),GAMA(NNNW)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /SCF1/UCF(NNNW)
     :      /SCF3/SCNSTY(NNNW),METHOD(NNNW)
     :      /SYMA/PIATJP,PIASPA
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
*
*   Write out the orbital properties
*
      WRITE (24,301)
      DO 1 I = 1,NW
         WRITE (24,302) NP(I),NH(I),E(I),PZ(I),
     :                  GAMA(I),PF(2,I),QF(2,I),SCNSTY(I),MF(I)
    1 CONTINUE
*
      WRITE (24,303)
      DO 2 I = 1,NW
         WA = RINT (I,I,-1)
         WB = RINT (I,I, 1)
         WC = RINT (I,I, 2)
         WD = RINT (I,I, 4)
         WE = 0.d0
         If (NH(I) /= 's ' .AND. NH(i) /= 'p-') then
            WE = RINT(I,I,-3)
         END IF
         WRITE (24,304) NP(I),NH(I),WE,WA,WB,WC,WD, UCF(I)
    2 CONTINUE
*
      IF (NCMIN .NE. 0) THEN
         mode = 0
         CALL ENGOUT (EVAL,IATJPO,IASPAR,ICCMIN,NCMIN,MODE)
         CALL CSFWGT (.FALSE.)
      ENDIF
*
      CLOSE (24)
*
      RETURN
*
  301 FORMAT (/'Radial wavefunction summary:'
     :       //67X, 'Self'
     :       /'Subshell',6X,'e',13X,'p0',5X,
     :         'gamma',5X,'P(2)',7X,'Q(2)',3X,'Consistency',' MTP'
     :        /)
  302 FORMAT (1X,I2,A2,1X,1P,D17.10,1P,D11.3,0P,F6.2,1P,3(D11.3),I5)
  303 FORMAT (/18X,'-3',14X,'-1',29X,'2',14x,'4',5X,'Generalised'
     :        /'Subshell',4X,'<  r  >',8X,'<  r  >',8X,'<  r  >',8X,
     :         '<  r  >',8X,'<  r  >',6X,'occupation'/)
  304 FORMAT (1X,1I2,1A2,1X,1P,6D15.5)
*
      END
