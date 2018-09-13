************************************************************************
*                                                                      *
      SUBROUTINE PRINTALS (INUM_II,INUM_FF,ASFA,ASFB,I,J,OMEGA,FACTOR)
*                                                                      *
*   This  routine  prints the basic oscillator strength  information   *
*   for transitions between level I and level J in LSJ coupling.       *
*                                                                      *
*   Written by G. Gaigalas,                                            *
*   NIST                                                  May 2011     *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      LOGICAL LTC
      CHARACTER*1 Lev_par_1,Lev_par_2,IM
      CHARACTER*4 Lev_J_1,Lev_J_2
      CHARACTER*64 string_CSF1,string_CSF2
      DIMENSION DBLFAC(10)
*
      POINTER(PIATJPII,IATJPOII(*)),(PIASPAII,IASPARII(*))
      POINTER(PIATJPFF,IATJPOFF(*)),(PIASPAFF,IASPARFF(*))
*
      COMMON/SYMAII/PIATJPII,PIASPAII
      COMMON/SYMAFF/PIATJPFF,PIASPAFF
*
      COMMON/CUTO/CUTOFF
     :      /DEF9/CVAC,PI
     :      /DEF10/AUCM,AUEV,CCMPS,FASI,FBSI
     :      /OSC2/LK,KK
     :      /OSC7/LTC(10)
     :      /DEF1/EMN,IONCTY,NELEC,Z
*
      POINTER(PNELev_POS_1,Lev_POS_1(*))
      POINTER(PNELev_J_1,Lev_J_1(*))
      POINTER(PNELev_Par_1,Lev_Par_1(*))
      POINTER(PNELev_ENER_1,RLev_ENER_1(*))
      POINTER(PNEstring_CSF1,string_CSF1(*))
*
      COMMON/JJ2LSJ1/ NVECTOTI,IOPEN_STATUS1,PNELev_POS_1,PNELev_J_1,
     :                PNELev_Par_1,PNELev_ENER_1,PNEstring_CSF1
*
      POINTER(PNELev_POS_2,Lev_POS_2(*))
      POINTER(PNELev_J_2,Lev_J_2(*))
      POINTER(PNELev_Par_2,Lev_Par_2(*))
      POINTER(PNELev_ENER_2,RLev_ENER_2(*))
      POINTER(PNEstring_CSF2,string_CSF2(*))
*
      COMMON/JJ2LSJ2/ NVECTOTF,IOPEN_STATUS2,PNELev_POS_2,PNELev_J_2,
     :                 PNELev_Par_2,PNELev_ENER_2,PNEstring_CSF2
*
*  *****  DBLFAC(I) = (2I+1)!!
*
      DATA DBLFAC/      3.0000000000D 00,1.5000000000D 01
     :,1.0500000000D 02,9.4500000000D 02,1.0395000000D 04
     :,1.3513500000D 05,2.0270250000D 06,3.4459425000D 07
     :,6.5472907500D 08,1.3749310575D 10
     :/
*
      D    = DABS(RLev_ENER_2(INUM_FF) - RLev_ENER_1(INUM_II))
      IF(DABS(D - DABS(OMEGA))  .GT. 0.00000001) THEN
       print*, '                                                      '
       print*, ' INCORRECT INPUT, PROGRAM STOP!!!                     '
       print*, '                                                      '
       print*, ' lsj.lbl files from the jj2lsj runs are inconsistent  ' 
       print*, ' with the files used in the bioscl calculation.       '
       print*, ' The bioscl calculation was for rci wave functions    '
       print*, ' but the jj2lsj run was for rscf wave functions or    '
       print*, ' vice versa.                                          '
       print*, '                                                      '
       print*, ' If the bioscl calculation was for RCI wave functions '
       print*, ' run jj2lsj to produce labels for RCI wave functions. '
       print*, ' Then rerun bioscl.                                   '
       print*, '                                                      '
       print*, ' If the bioscl calculation was for rscf wave functions'
       print*, ' run jj2lsj to produce labels for rscf wave functions.'
       print*, ' Then rerun bioscl.                                   '
       STOP
      END IF
CGG      DD    = D*AUCM*(1.0+1.0/EMN) 
      DD    = D*AUCM 
      ANGS  =  1.0D 08 / DD 
      ANGSA = ANGS
      IF(ANGS .GT. 2000.D0) THEN
         SIGMA = (1.D8/ANGS)**2
         ANGSA = ANGS/(1.D0 + 8342.13D-8 + 206030./(130.D+8 - SIGMA)
     :                 +15997./(38.9D+8 - SIGMA))
      END IF
*
*   Evaluate statistical factors and constants
*
      STFAC = IATJPOII(I)
      DLL1 = DBLE (LK+LK+1)
      STFAC = STFAC/DLL1
*
      OMC = OMEGA/CVAC
      FAAU = 2.0D 00*OMC*STFAC/DBLE (IATJPOFF(J))
      FOSC = CVAC*STFAC/OMC
      ENG = OMEGA*FACTOR
      IF (LTC(1)) ENG = FACTOR/OMEGA
*
      GGFACTOR = CVAC**(2*LK-2)*DBLE(LK)*DBLFAC(LK)**2/
     :           ((2.0D 00 * DBLE(LK)+1.0D 00)*
     :           (DBLE(LK)+1.0D 00)*ABS(OMEGA)**(2*LK-1))
*
*   Calculate Einstein A and B coefficients and oscillator strengths
*
      IF (KK .EQ. 0) THEN
*
*   Electric multipoles
*
*   In atomic units
*
         ACSQ = ASFA**2
         ABSQ = ASFB**2
         AC = ACSQ*FAAU
         AB = ABSQ*FAAU
         OSCC = ACSQ*FOSC
         OSCB = ABSQ*FOSC
         SA=OSCC*GGFACTOR
         SB=OSCB*GGFACTOR
*
*   Convert to SI units if option 5 not set
*
         IF (.NOT. LTC(7)) THEN
            AC = AC*FASI
            AB = AB*FASI
         ENDIF
*
*   Print information if both AC and AB are greater than CUTOFF
*
         IF((ABS(AC).GE.CUTOFF) .AND. (ABS(AB).GE.CUTOFF)) THEN
            IM = 'E'
            IF (-1.d+9 < ENG .and. ENG .LT. -1.0d-2) THEN
               WRITE (32,'(A1)')  ""
               WRITE (32,'(A1)')  ""
               IF(DABS(RLev_ENER_1(INUM_II)) .LT. 9999.5) THEN
                  WRITE(32,'(I4,F14.8,2X,A)')
     :            IATJPOII(I)-1,RLev_ENER_1(INUM_II), 
     :            trim(string_CSF1(INUM_II))
               ELSE IF(DABS(RLev_ENER_1(INUM_II)) .LT. 99999.5) THEN
                   WRITE(32,'(I4,F14.7,2X,A)')
     :             IATJPOII(I)-1,RLev_ENER_1(INUM_II), 
     :             trim(string_CSF1(INUM_II))
               ELSE
                   WRITE(32,'(I4,F14.6,2X,A)')
     :             IATJPOII(I)-1,RLev_ENER_1(INUM_II), 
     :             trim(string_CSF1(INUM_II))
               END IF
               IF(DABS(RLev_ENER_2(INUM_FF)) .LT. 9999.5) THEN
                  WRITE(32,'(I4,F14.8,2X,A)')
     :            IATJPOFF(J)-1,RLev_ENER_2(INUM_FF), 
     :            trim(string_CSF2(INUM_FF))
               ELSE IF(DABS(RLev_ENER_2(INUM_FF)) .LT. 99999.5) THEN
                  WRITE(32,'(I4,F14.7,2X,A)')
     :            IATJPOFF(J)-1,RLev_ENER_2(INUM_FF), 
     :            trim(string_CSF2(INUM_FF))
               ELSE
                  WRITE(32,'(I4,F14.6,2X,A)')
     :            IATJPOFF(J)-1,RLev_ENER_2(INUM_FF), 
     :            trim(string_CSF2(INUM_FF))
               END IF
               WRITE (32,38) DD,ANGS,ANGSA,
     :             IM,LK,-SB,-OSCB,-AB,
     :             DABS(DABS(SB)-DABS(SA))/MAX(DABS(SB),DABS(SA))
CGG     :             IM,LK,-SB,-OSCB,-AB
               WRITE (32,39)  -SA,-OSCC,-AC
            ELSE IF ( 1.d+9 .GT. ENG .AND. ENG .GT. 1.D-2) THEN
               WRITE (32,'(A1)')  ""
               WRITE (32,'(A1)')  ""
               IF(DABS(RLev_ENER_2(INUM_FF)) .LT. 9999.5) THEN
                  WRITE(32,'(I4,F14.8,2X,A)')
     :            IATJPOFF(J)-1,RLev_ENER_2(INUM_FF), 
     :            trim(string_CSF2(INUM_FF))
               ELSE IF(DABS(RLev_ENER_2(INUM_FF)) .LT. 99999.5) THEN
                  WRITE(32,'(I4,F14.7,2X,A)')
     :            IATJPOFF(J)-1,RLev_ENER_2(INUM_FF), 
     :            trim(string_CSF2(INUM_FF))
               ELSE
                  WRITE(32,'(I4,F14.6,2X,A)')
     :            IATJPOFF(J)-1,RLev_ENER_2(INUM_FF), 
     :            trim(string_CSF2(INUM_FF))
               END IF
               IF(DABS(RLev_ENER_1(INUM_II)) .LT. 9999.5) THEN
                  WRITE(32,'(I4,F14.8,2X,A)')
     :            IATJPOII(I)-1,RLev_ENER_1(INUM_II), 
     :            trim(string_CSF1(INUM_II))
               ELSE IF(DABS(RLev_ENER_1(INUM_II)) .LT. 99999.5) THEN
                   WRITE(32,'(I4,F14.7,2X,A)')
     :             IATJPOII(I)-1,RLev_ENER_1(INUM_II), 
     :             trim(string_CSF1(INUM_II))
               ELSE
                   WRITE(32,'(I4,F14.6,2X,A)')
     :             IATJPOII(I)-1,RLev_ENER_1(INUM_II), 
     :             trim(string_CSF1(INUM_II))
               END IF
               WRITE (32,38) DD,ANGS,ANGSA,
     :             IM,LK,SB,OSCB,AB*IATJPOFF(J)/IATJPOII(I),
     :             DABS(DABS(SB)-DABS(SA))/MAX(DABS(SB),DABS(SA))
CGG     :             IM,LK,SB,OSCB,AB*IATJPOFF(J)/IATJPOII(I)
               WRITE (32,39)  SA,OSCC,AC*IATJPOFF(J)/IATJPOII(I)
            END IF
         ENDIF
      ELSE
*
*   Magnetic multipoles
*
*   In atomic units
*
         AMS = ASFA**2
         AM = AMS*FAAU
         OSCM = AMS*FOSC
         SA=OSCM*GGFACTOR*CVAC*CVAC*4
*
*   Convert to SI units if option 5 not set
*
         IF (.NOT. LTC(7)) AM = AM*FASI
*
*   Print information if AM is greater than CUTOFF
*
         IF (ABS(AM) .GE. CUTOFF) THEN
            IM = 'M'
            IF (-1.d+9 < ENG .and. ENG .LT. -1.0d-2) THEN
               WRITE (32,'(A1)')  ""
               WRITE (32,'(A1)')  ""
               IF(DABS(RLev_ENER_1(INUM_II)) .LT. 9999.5) THEN
                  WRITE(32,'(I4,F14.8,2X,A)')
     :            IATJPOII(I)-1,RLev_ENER_1(INUM_II), 
     :            trim(string_CSF1(INUM_II))
               ELSE IF(DABS(RLev_ENER_1(INUM_II)) .LT. 99999.5) THEN
                   WRITE(32,'(I4,F14.7,2X,A)')
     :             IATJPOII(I)-1,RLev_ENER_1(INUM_II), 
     :             trim(string_CSF1(INUM_II))
               ELSE
                   WRITE(32,'(I4,F14.6,2X,A)')
     :             IATJPOII(I)-1,RLev_ENER_1(INUM_II), 
     :             trim(string_CSF1(INUM_II))
               END IF
               IF(DABS(RLev_ENER_2(INUM_FF)) .LT. 9999.5) THEN
                  WRITE(32,'(I4,F14.8,2X,A)')
     :            IATJPOFF(J)-1,RLev_ENER_2(INUM_FF), 
     :            trim(string_CSF2(INUM_FF))
               ELSE IF(DABS(RLev_ENER_2(INUM_FF)) .LT. 99999.5) THEN
                  WRITE(32,'(I4,F14.7,2X,A)')
     :            IATJPOFF(J)-1,RLev_ENER_2(INUM_FF), 
     :            trim(string_CSF2(INUM_FF))
               ELSE
                  WRITE(32,'(I4,F14.6,2X,A)')
     :            IATJPOFF(J)-1,RLev_ENER_2(INUM_FF), 
     :            trim(string_CSF2(INUM_FF))
               END IF
               WRITE (32,40) DD,ANGS,ANGSA,IM,LK,-SA,-OSCM,-AM
CGG     :             IM,LK,-SA,-OSCM,-AM
            ELSE IF ( 1.d+9 .GT. ENG .AND. ENG .GT. 1.D-2) THEN
               WRITE (32,'(A1)')  ""
               WRITE (32,'(A1)')  ""
               IF(DABS(RLev_ENER_2(INUM_FF)) .LT. 9999.5) THEN
                  WRITE(32,'(I4,F14.8,2X,A)')
     :            IATJPOFF(J)-1,RLev_ENER_2(INUM_FF), 
     :            trim(string_CSF2(INUM_FF))
               ELSE IF(DABS(RLev_ENER_2(INUM_FF)) .LT. 99999.5) THEN
                  WRITE(32,'(I4,F14.7,2X,A)')
     :            IATJPOFF(J)-1,RLev_ENER_2(INUM_FF), 
     :            trim(string_CSF2(INUM_FF))
               ELSE
                  WRITE(32,'(I4,F14.6,2X,A)')
     :            IATJPOFF(J)-1,RLev_ENER_2(INUM_FF), 
     :            trim(string_CSF2(INUM_FF))
               END IF
               IF(DABS(RLev_ENER_1(INUM_II)) .LT. 9999.5) THEN
                  WRITE(32,'(I4,F14.8,2X,A)')
     :            IATJPOII(I)-1,RLev_ENER_1(INUM_II), 
     :            trim(string_CSF1(INUM_II))
               ELSE IF(DABS(RLev_ENER_1(INUM_II)) .LT. 99999.5) THEN
                   WRITE(32,'(I4,F14.7,2X,A)')
     :             IATJPOII(I)-1,RLev_ENER_1(INUM_II), 
     :             trim(string_CSF1(INUM_II))
               ELSE
                   WRITE(32,'(I4,F14.6,2X,A)')
     :             IATJPOII(I)-1,RLev_ENER_1(INUM_II), 
     :             trim(string_CSF1(INUM_II))
               END IF
               WRITE (32,40) DD,ANGS,ANGSA,
     :             IM,LK,SA,OSCM,AM*IATJPOFF(J)/IATJPOII(I)
CGG     :             IM,LK,SA,OSCM,AM*IATJPOFF(J)/IATJPOII(I)
            END IF 
         ENDIF
      ENDIF
      RETURN
   38 FORMAT(F11.2,' CM-1',2X,F11.2,' ANGS(VAC)',2X,F11.2,' ANGS(AIR)'/
     :    1X,A1,I1,2X,'S = ',1PD12.5,3X,'GF = ',D12.5,3X,'AKI = ',D12.5,
     :    3X,'dT = ',0PF8.5)
CGG     :    1X,A1,I1,2X,'S = ',1PD12.5,3X,'GF = ',D12.5,3X,'AKI = ',D12.5)
   39 FORMAT(9X,1PD12.5,8X,D12.5,9X,D12.5)
   40 FORMAT(F11.2,' CM-1',2X,F11.2,' ANGS(VAC)',2X,F11.2,' ANGS(AIR)'/
     :    1X,A1,I1,2X,'S = ',1PD12.5,3X,'GF = ',D12.5,3X,'AKI = ',D12.5)
      END
