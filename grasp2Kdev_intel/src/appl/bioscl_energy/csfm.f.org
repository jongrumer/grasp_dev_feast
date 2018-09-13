************************************************************************
*                                                                      *
      SUBROUTINE CSFM (ASFA,ASFB,LEV1,LEV2)
*                                                                      *
*   This routine calculates  the CSF Coulomb, Babuskin, and magnetic   *
*   matrix elements for  a transition  between  levels  separated by   *
*   energy OMEGA.                                                      *
*                                                                      *
*   Modified for different initial and final state orbitals            *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)

      include 'parameters.def'
Cff      PARAMETER (NNNW = 120)
Cff      PARAMETER (KEY = 121)
      PARAMETER (KEY = KEYORB)
      PARAMETER (NCA = 65536)

Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL LDBPR
      CHARACTER*2 NH
      CHARACTER*4 JLBL,LABJ,LABP,JLABI,JLABJ,JPARI,JPARJ                       
*
      POINTER (PNEVEC,EVEC(1))
      POINTER (PNEVECII,EVECII(1))
      POINTER (PNEVECFF,EVECFF(1))

      POINTER(PIATJPII,IATJPOII(1)),(PIASPAII,IASPARII(1))
      POINTER(PIATJPFF,IATJPOFF(1)),(PIASPAFF,IASPARFF(1))                      

      POINTER (PNTJJA,JJA(1)),(PNTJJB,JJB(1)),
     :        (PNTHB1,HB1(1)),(PNTHC1,HB2(1)),
     :        (PNTHB2,HC1(1)),(PNTHC2,HC2(1)),
     :        (PNTHM1,HM1(1)),(PNTHM2,HM2(1))
      POINTER (PISLDR,ISLDR(1)),(PXSLDR,XSLDR(1))
      POINTER (PISLDR1,ISLDR1(1))
      POINTER (PNTLAB,LAB(1)),(PNNPTR,NPTR(1))
*
      COMMON/SYMAII/PIATJPII,PIASPAII
     :      /PRNTII/NVECII,PNIVECII,NVECMXII
 
      COMMON/SYMAFF/PIATJPFF,PIASPAFF
     :      /PRNTFF/NVECFF,PNIVECFF,NVECMXFF                                    
*
      COMMON/DEBUGR/LDBPR(30),CUTOFF
     :      /EIGVEC/PNEVEC
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
     :      /OSC1/PNTJJA,PNTJJB,PNTHB1,PNTHB2,
     :            PNTHC1,PNTHC2,PNTHM1,PNTHM2,NSDIM
     :      /OSC2/LK,KK
     :      /OSC3/PXSLDR,PISLDR,PISLDR1,NTDIM
     :      /OSC5/NINT,PNTLAB,PNNPTR,NINTEG
     :      /JLABL/JLBL(32),LABJ(32),LABP(2)
     :      /SYMA/PIATJP,PIASPA                                                 

      COMMON/EIGVECII/PNEVECII
     :      /ORB2II/NCFII,NWII

      COMMON/EIGVECFF/PNEVECFF
     :      /ORB2FF/NCFFF,NWFF

       POINTER (PNIVECII,IVECII(1))
       POINTER (PNIVECFF,IVECFF(1))                                             

*
      EPS = 1.0D-10
*
      ASFA = 0.0D 00
      ASFB = 0.0D 00
      NCOUNT = 0
      NFI = 0
*
*Rasa -- start
*
*   J/pi labels for levels
*
      JLABI = LABJ(IATJPOII(LEV1))
      JLABJ = LABJ(IATJPOFF(LEV2))
      IPAR = (IASPARII(LEV1)+3)/2
      JPAR = (IASPARFF(LEV2)+3)/2
      JPARI = LABP(IPAR)
      JPARJ = LABP(JPAR)                                                        
     
      IF (LDBPR(18)) THEN
          WRITE (*,302) IVECFF(LEV2), JLABJ, JPARJ, IVECII(LEV1), 
     * JLABI, JPARI
          WRITE (*,301)
      ENDIF
* Rasa -- end
      DO M = 1,NINTEG
        JA = LAB(M)
        JB = MOD (JA,KEY)
        JA = JA/KEY
        NCLR = NPTR(M)+1
        NCUP = NPTR(M+1)
        IF (NCLR .GT. NCUP) GOTO 4
        IF (MOD(NKL(JA)+NKL(JB)+KK+LK,2) .NE. 0) GOTO 4
        CALL SPME (JA,JB,HCOUL1,HBAB1,HMAG1)
*
        DO J = NCLR,NCUP
*
          IA = ISLDR(J)
          IB = ISLDR1(J)    
c         IB = MOD (IA,NCA)
c         IA = IA/NCA
*
*  Check on consistency
*
          IF (IA.GT.NCFII) THEN
            GOTO 131
          ELSEIF (IB.GT.NCFFF) THEN
            GOTO 131
          ENDIF
*
*  Observ that IA and IB refer to the merged list
*  whereas EVECII and EVECFF refers to the initial and
*  final state lists. Therefore we have IB -> IB-NCFII
*
          COUVX =EVECII(IA+(LEV1-1)*NCFII)*
     :           EVECFF(IB+(LEV2-1)*NCFFF)
c     :           EVECFF(IB-NCFII+(LEV2-1)*NCFFF)
          COEFF = XSLDR(J)
          IF (ABS(COUVX) .GT. EPS) THEN
CRasa       IF (LDBPR(18)) WRITE (99,300) NP(JA),NH(JA),
CRasa:      NP(JB),NH(JB),IA,IB,COEFF
            IF (KK .EQ. 0) THEN
              ASFA = ASFA+HCOUL1*COEFF*COUVX
              ASFB = ASFB+HBAB1*COEFF*COUVX
*Rasa -- start
              IF (LDBPR(18) .and. (dabs(HCOUL1*COEFF*COUVX)
     *          .gt. cutoff)) THEN
                  WRITE (*,300) IA, EVECII(IA+(LEV1-1)*NCFII),
     *                IB-NCFII, EVECFF(IB-NCFII+(LEV2-1)*NCFFF),
     *                NP(JA),NH(JA),NP(JB),NH(JB),COEFF,HCOUL1,'C',
     *                HCOUL1*COEFF*COUVX
                  WRITE(*,303) HBAB1,'B',HBAB1*COEFF*COUVX
              ENDIF
*Rasa -- end
            ELSE
              contr=HMAG1*COEFF*COUVX
              ASFA = ASFA+HMAG1*COEFF*COUVX
*Rasa -- start
              IF (LDBPR(18) .and. (dabs(contr) .gt. CUTOFF)) THEN
                  WRITE (*,300) IA, EVECII(IA+(LEV1-1)*NCFII),
     *                IB-NCFII, EVECFF(IB-NCFII+(LEV2-1)*NCFFF),
     *                NP(JA),NH(JA),NP(JB),NH(JB),COEFF,HMAG1,' ',
     *                contr
              ENDIF
*Rasa -- end
            ENDIF
          ENDIF
 131    CONTINUE
        ENDDO
    4 ENDDO
      Write(*,304)
      IF (KK .EQ. 0) THEN
          WRITE(*,305) 'C',ASFA
          WRITE(*,305) 'B',ASFB
      else
          WRITE(*,305) ' ',ASFA
      endif
*
      RETURN
*
  300 FORMAT (I5,2X,D11.5,I5,2X,D11.5, 2(1X,1I2,1A2),1P,D15.6,1X,D15.6,
     *1X,1A,1X,D15.6)
  301 FORMAT (4X,'i',3X,'Coeff.(i)',5X,'f',2X,'Coeff.(f)'2X,
     *'Orb(i)',1X,'Orb(f)',2X,'MCT coeff.',4X,'Radial Factor',
     *5X,'Contribution')
  302 FORMAT ('Upper level ',I3,A4,A4,10X,' Lower level',I3,A4,A4)
          WRITE (*,302) IVECFF(LEV2), JLABJ, JPARJ, LEV1, JLABI, JPARI
  303 FORMAT (62X,1P,D15.6,1X,1A,1X,D15.6)
  304 FORMAT(65X,'Total:')
  305 FORMAT (62X,1P,16X,1A,1X,D15.6)
*
      END
