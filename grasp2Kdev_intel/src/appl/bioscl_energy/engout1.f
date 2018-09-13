************************************************************************
*                                                                      *
      SUBROUTINE ENGOUT1 (EAV,E,JTOT,IPAR,ILEV,NN,MODE,K)
*                                                                      *
*   This  subroutine prints  energy levels, splittings, and energies   *
*   relative to the lowest in  Hartrees, Kaysers, and  eV, using the   *
*   reduced mass corrected value for the Rydberg. If MODE is 0, only   *
*   the eigenenergies are printed. If  MODE  is 1, the eigenenergies   *
*   and separations are printed. If MODE is 2, the eigenenergies and   *
*   energies relative to level 1 are printed. If MODE is 3, the eig-   *
*   enenergies,  separations,  and energies relative to level  1 are   *
*   printed.                                                           *
*                                          Last updated: 15 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      CHARACTER*4 JLBL,LABJ,LABP
*
      DIMENSION E(NN),JTOT(NN),IPAR(NN),ILEV(NN)
*
      COMMON/DEF10/AUCM,AUEV,CCMS,FASI,FBSI
     :      /JLABL/JLBL(32),LABJ(32),LABP(2)
*
*   Always print the eigenenergies
*
      IF (K.EQ.1) WRITE (24,299)
      IF (K.EQ.2) WRITE (24,300)
      WRITE (24,301)
      DO 1 J = 1,NN
         I = ILEV(J)
         EAU = E(J)+EAV
         ECM = EAU*AUCM
         EEV = EAU*AUEV
         IP = (IPAR(J)+3)/2
         WRITE (24,302) I,LABJ(JTOT(J)),LABP(IP),EAU,ECM
    1 CONTINUE
*****************************************************************
* Commented out on 98-10-16
*      IF (NN .GT. 1) THEN
*
*   Energy separations
*
*         IF ((MODE .EQ. 1) .OR. (MODE .EQ. 3)) THEN
*            WRITE (24,303)
*            WRITE (24,301)
*            DO 2 J = 2,NN
*               I = ILEV(J)
*               EAU = E(J)-E(J-1)
*               ECM = EAU*AUCM
*               EEV = EAU*AUEV
*               IP = (IPAR(J)+3)/2
*               WRITE (24,302) I,LABJ(JTOT(J)),LABP(IP),EAU,ECM
*    2       CONTINUE
*         ENDIF
**
**   Energies relative to level 1
**
*         IF ((MODE .EQ. 2) .OR. (MODE .EQ. 3)) THEN
*            WRITE (24,304)
*            WRITE (24,301)
*            DO 3 J = 2,NN
*               I = ILEV(J)
*               EAU = E(J)-E(1)
*               ECM = EAU*AUCM
*               EEV = EAU*AUEV
*               IP = (IPAR(J)+3)/2
*               WRITE (24,302) I,LABJ(JTOT(J)),LABP(IP),EAU,ECM
*    3       CONTINUE
*         ENDIF
*
*      ENDIF
*****************************************************************
*
      RETURN
*
  299 FORMAT ('Eigenenergies for the initial state list')
  300 FORMAT ('Eigenenergies for the final state list')
  301 FORMAT ('Level  J Parity',10X,'Hartrees',18X,'Kaysers')
  302 FORMAT (1I3,4X,2A4,1P,2D25.15)
  303 FORMAT ('Energy of each level relative to immediately lower',
     :          ' level:')
  304 FORMAT ('Energy of each level relative to lowest level:')
*
      END
