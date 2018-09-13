************************************************************************
*                                                                      *
      SUBROUTINE ENGOUTH (EAV,E,JTOT,IPAR,ILEV,NN,MODE)
*                                                                      *
*   This  subroutine prints  energy levels, splittings, and energies   *
*   relative to the lowest in  Hartrees, Kaysers, and  eV, using the   *
*   reduced mass corrected value for the Rydberg. If MODE is 0, only   *
*   the eigenenergies are printed. If  MODE  is 1, the eigenenergies   *
*   and separations are printed. If MODE is 2, the eigenenergies and   *
*   energies relative to level 1 are printed. If MODE is 3, the eig-   *
*   enenergies,  separations,  and energies relative to level  1 are   *
*   printed.                                                           *
*                                          Last updated: 27 July 2006  *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      CHARACTER*4 JLBL,LABJ,LABP
*
      DIMENSION E(NN),JTOT(NN),IPAR(NN),ILEV(NN)
*
      COMMON/DEF10/AUCM,AUEV,CCMS,FASI,FBSI
     :      /JLABL/JLBL(32),LABJ(32),LABP(2)
*
*   Always print the eigenenergies
*
      WRITE (29,300)
      WRITE (29,301)
      DO J = 1,NN
         I = ILEV(J)
         EAU = E(J)+EAV
         ECM = EAU*AUCM
         EEV = EAU*AUEV
         IP = (IPAR(J)+3)/2
         WRITE (29,302) I,LABJ(JTOT(J)),LABP(IP),EAU,ECM,EEV
      ENDDO
*
      IF (NN .GT. 1) THEN
*
*   Energy separations
*
         IF ((MODE .EQ. 1) .OR. (MODE .EQ. 3)) THEN
            WRITE (29,303)
            WRITE (29,301)
            DO J = 2,NN
               I = ILEV(J)
               EAU = E(J)-E(J-1)
               ECM = EAU*AUCM
               EEV = EAU*AUEV
               IP = (IPAR(J)+3)/2
               WRITE (29,302) I,LABJ(JTOT(J)),LABP(IP),EAU,ECM,EEV
            ENDDO
         ENDIF
*
*   Energies relative to level 1
*
         IF ((MODE .EQ. 2) .OR. (MODE .EQ. 3)) THEN
            WRITE (29,304)
            WRITE (29,301)
            DO J = 2,NN
               I = ILEV(J)
               EAU = E(J)-E(1)
               ECM = EAU*AUCM
               EEV = EAU*AUEV
               IP = (IPAR(J)+3)/2
               WRITE (29,302) I,LABJ(JTOT(J)),LABP(IP),EAU,ECM,EEV
            ENDDO
         ENDIF
*
      ENDIF
*
      RETURN
*
  300 FORMAT (/'Eigenenergies:')
  301 FORMAT (/'Level  J Parity',7X,'Hartrees',14X,'Kaysers',
     :         16X,'eV'/)
  302 FORMAT (1I3,2X,2A4,1P,3D22.14)
  303 FORMAT (/'Energy of each level relative to immediately lower',
     :          ' level:')
  304 FORMAT (/'Energy of each level relative to lowest level:')
  305 FORMAT (3X,'J',8X,'Hartrees')
  306 FORMAT (A4,2X,D22.14)
*
      END
