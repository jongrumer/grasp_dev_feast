************************************************************************
*                                                                      *
      SUBROUTINE PRINTA (ASFA,ASFB,I,J,OMEGA,FACTOR,LINES,LSAME)
*                                                                      *
*   This  routine  prints the basic oscillator strength  information   *
*   for transitions between level I and level J.                       *
*                                                                      *
*                                         Last revision: 28 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      LOGICAL LTC,LSAME
      CHARACTER*4 JLBL,LABJ,LABP,JLABI,JLABJ,JPARI,JPARJ
      CHARACTER*2 F1,F2
CGG
      DIMENSION DBLFAC(10)
CGG-end
*
      POINTER (PIATJP,IATJPO(*)),(PIASPA,IASPAR(*))
      POINTER(PIATJPII,IATJPOII(*)),(PIASPAII,IASPARII(*))
      POINTER(PIATJPFF,IATJPOFF(*)),(PIASPAFF,IASPARFF(*))
      POINTER (PNTOTB,TOTB(*)),(PNTOTC,TOTC(*))
c
cbieron  numbering of levels in OSCL
c
       POINTER (PNIVEC,IVEC(*))
       POINTER (PNIVECII,IVECII(*))
       POINTER (PNIVECFF,IVECFF(*))
*
      COMMON/CUTO/CUTOFF
     :      /DEF2/C
     :      /DEF9/CVAC,PI
     :      /DEF10/AUCM,AUEV,CCMPS,FASI,FBSI
     :      /JLABL/JLBL(32),LABJ(32),LABP(2)
     :      /SYMA/PIATJP,PIASPA
     :      /OSC2/LK,KK
     :      /OSC4/PNTOTC,PNTOTB
     :      /OSC7/LTC(10)

      COMMON/SYMAII/PIATJPII,PIASPAII
     :      /PRNTII/NVECII,PNIVECII,NVECMXII

      COMMON/SYMAFF/PIATJPFF,PIASPAFF
     :      /PRNTFF/NVECFF,PNIVECFF,NVECMXFF
c
cbieron  numbering of levels in OSCL
c
     :      /PRNT/NVEC,PNIVEC,NVECMX
CGG
*
*  *****  DBLFAC(I) = (2I+1)!!
*
      DATA DBLFAC/      3.0000000000D 00,1.5000000000D 01
     :,1.0500000000D 02,9.4500000000D 02,1.0395000000D 04
     :,1.3513500000D 05,2.0270250000D 06,3.4459425000D 07
     :,6.5472907500D 08,1.3749310575D 10
     :/
CGG-end
*
*   Evaluate statistical factors and constants
*
      STFAC = IATJPOII(I)
      DLL1 = DBLE (LK+LK+1)
      STFAC = STFAC/DLL1
*
      OMC = OMEGA/CVAC

      print*, 'OMEGA = ', OMEGA

c     OMC = OMEGA/C
      FAAU = 2.0D 00*OMC*STFAC/DBLE (IATJPOFF(J))
      FOSC = CVAC*STFAC/OMC
c     FOSC = C*STFAC/OMC
      FBAU = PI*FOSC/OMEGA
      ENG = OMEGA*FACTOR
      IF (LTC(1)) ENG = FACTOR/OMEGA
*
*   J/pi labels for levels
*
      JLABI = LABJ(IATJPOII(I))
      JLABJ = LABJ(IATJPOFF(J))
      IPAR = (IASPARII(I)+3)/2
      JPAR = (IASPARFF(J)+3)/2
      JPARI = LABP(IPAR)
      JPARJ = LABP(JPAR)
CGG
      GGFACTOR = CVAC**(2*LK-2)*DBLE(LK)*DBLFAC(LK)**2/
     :((2.0D 00 * DBLE(LK)+1.0D 00)*
     : (DBLE(LK)+1.0D 00)*ABS(OMEGA)**(2*LK-1))
CGG-end
*
*   Calculate Einstein A and B coefficients and oscillator strengths
*
      IF (KK .EQ. 0) THEN
*
*   Electric multipoles
*
*   In atomic units
*
         IF (ASFA.LT.0) THEN
            ISIGNA = -1
         ELSE
            ISIGNA = 1
         ENDIF
         IF (ASFB.LT.0) THEN
            ISIGNB = -1
         ELSE
            ISIGNB = 1
         ENDIF

         ACSQ = ASFA**2
         ABSQ = ASFB**2
         AC = ACSQ*FAAU
         AB = ABSQ*FAAU
         BC = ACSQ*FBAU
         BB = ABSQ*FBAU
         OSCC = ACSQ*FOSC
         OSCB = ABSQ*FOSC
CGG
CGG         SA=OSCC*1.5/ABS(OMEGA)
CGG         SB=OSCB*1.5/ABS(OMEGA)
         SA=OSCC*GGFACTOR
         SB=OSCB*GGFACTOR
CGG-end
*
*   Convert to SI units if option 5 not set
*
         IF (.NOT. LTC(7)) THEN
            AC = AC*FASI
            AB = AB*FASI
            BC = BC*FBSI
            BB = BB*FBSI
         ENDIF
*
*   Accumulate total of A coefficients
*
         TOTC(J) = TOTC(J)+AC
         TOTB(J) = TOTB(J)+AB
*
*   Print information if both AC and AB are greater than CUTOFF
*
         IF ((ABS (AC) .GE. CUTOFF) .AND. (ABS (AB) .GE. CUTOFF))
     :      THEN
c
cbieron  numbering of levels in OSCL
c
c            WRITE (24,300) J,JLABJ,JPARJ,I,JLABI,JPARI,ENG,
!xhh            WRITE (24,300) IVECFF(J),JLABJ,JPARJ,I,JLABI,JPARI,ENG,

CPer 2007-03-24 Write transitions in the following way:
C       Energy of the upper level should be greater than the energy of
C       the lower level
C       Energy difference will thus always be positive
C       A in emission so that T = 1/A
C       All quantities will be positive
CFF     .. Furthermore
C       Transition energy is in units of kayzers
C       Print only if 0.01 < | ENG | < 1.0E+08

            IF (-1.d+9 < ENG .and. ENG .LT. -1.0d-2) THEN
               IF (LSAME) THEN
                  F1 = 'f '
                  F2 = 'f '
               ELSE
                  F1 = 'f2'
                  F2 = 'f1'
               END IF 
               WRITE (24,300) F1,IVECFF(J),JLABJ,JPARJ,F2,IVECII(I),
     :               JLABI,JPARI,-ENG,-AC,-OSCC,-SA
c     :              ,ASFA !(M-value)
C              WRITE (24,*) 'Relative sign',ISIGNA
               WRITE (24,301) -AB,-OSCB,-SB
c     :              ,ASFB !(M-value)
C              WRITE (24,*) 'Relative sign',ISIGNB
            ELSE IF ( 1.d+9 .GT. ENG .AND. ENG .GT. 1.D-2) THEN
               IF (LSAME) THEN
                  F1 = 'f '
                  F2 = 'f '
               ELSE
                  F1 = 'f1'
                  F2 = 'f2'
               END IF
               WRITE (24,300) F1,IVECII(I),JLABI,JPARI,F2,IVECFF(J),
     :               JLABJ,JPARJ,ENG,AC*IATJPOFF(J)/IATJPOII(I),
     :               OSCC,SA
c     :               ,ASFA !(M value)
C              WRITE (24,*) 'Relative sign',ISIGNA
               WRITE (24,301) AB*IATJPOFF(J)/IATJPOII(I),OSCB,SB
c     :               ,ASFB !(M value)
C              WRITE (24,*) 'Relative sign',ISIGNB
            END IF
            LINES = LINES+3
         ENDIF
*
      ELSE
*
*   Magnetic multipoles
*
*   In atomic units
*
         IF (ASFA.LT.0) THEN
            ISIGNA = -1
         ELSE
            ISIGNA = 1
         ENDIF
         AMS = ASFA**2
         AM = AMS*FAAU
         BM = AMS*FBAU
         OSCM = AMS*FOSC
CGG
CGG         SA=OSCM*1.5/ABS(OMEGA)
         SA=OSCM*GGFACTOR*CVAC*CVAC*4
CGG-end
*
*   Convert to SI units if option 5 not set
*
         IF (.NOT. LTC(7)) THEN
            AM = AM*FASI
            BM = AM*FBSI
         ENDIF
*
*   Accumulate total of A coefficients
*
         TOTC(J) = TOTC(J)+AM
         TOTB(J) = TOTB(J)+AM
*
*   Print information if AM is greater than CUTOFF
*
         IF (ABS (AM) .GE. CUTOFF) THEN
c
cbieron  numbering of levels in OSCL
c
c            WRITE (24,302) J,JLABJ,JPARJ,I,JLABI,JPARI,ENG,
            IF (-1.d+9 < ENG .and. ENG .LT. -1.0d-2) THEN
               IF (LSAME) THEN
                  F1 = 'f '
                  F2 = 'f '
               ELSE
                  F1 = 'f2'
                  F2 = 'f1'
               END IF 
               WRITE (24,302) F1,IVECFF(J),JLABJ,JPARJ,F2,IVECII(I),
     :               JLABI,JPARI,-ENG,-AM,-OSCM,-SA
C     :               ,ASFA (M-value)
C            WRITE (24,*) 'Relative sign',ISIGNA
            ELSE IF ( 1.d+9 .GT. ENG .AND. ENG .GT. 1.D-2) THEN
               IF (LSAME) THEN
                  F1 = 'f '
                  F2 = 'f '
               ELSE
                  F1 = 'f1'
                  F2 = 'f2'
               END IF
               WRITE (24,302) F1,IVECII(I),JLABI,JPARI,F2,IVECFF(J),
     :               JLABJ,JPARJ,ENG,AM*IATJPOFF(J)/IATJPOII(I),
     :               OSCM,SA
C     :              ,ASFA (M-value)
            END IF 
            LINES = LINES+2
         ENDIF
*
      ENDIF
*
      RETURN
*
  300 FORMAT(1X,A2,I3,1X,2A4,A2,I3,1X,2A4,0P,F11.2,' C',1P,
     :     3D13.5)
  301 FORMAT(40X,' B',1P,3D13.5)
  302 FORMAT(1X,A2,I3,1X,2A4,A2,I3,1X,2A4,0P,F11.2,' M',1P,
     :     3D13.5)
*
      END
