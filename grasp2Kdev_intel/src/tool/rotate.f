************************************************************************
*                                                                      *
      PROGRAM ROTATE                   
*                                                                      *
*   This program rotate orbitals with the same kappa                   *
*                                                                      *
*   Written by Per Jonsson,  Vanderbilt University 16 June 1996        *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = 600)
CGG      PARAMETER (NNNW = 120)

      INTEGER PNTRIQ
      LOGICAL GETYN, YES
      CHARACTER*2 NH
      CHARACTER*24 NAME
      CHARACTER*256 FILNAM
*
      POINTER (PNTRPF,PF(NNNP,1)),(PNTRQF,QF(NNNP,1))
*
      COMMON/GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /ORB1/E(NNNW),GAMA(NNNW)
     :      /ORB2/NCF,NW,PNTRIQ
      COMMON/ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
      COMMON/DEFAULT/NDEF
*
      WRITE(*,*) ' Welcome to the rotate program'
      WRITE(*,*) 
      WRITE(*,*) ' Input : name'
      WRITE(*,*) ' Output: name.rot'
      WRITE(*,*) 
      WRITE(*,*) ' Default settings?'
      YES = GETYN ()
      PRINT *
      IF (YES) THEN
         NDEF = 0
      ELSE
         NDEF = 1
      ENDIF

   10 PRINT *, 'Name of state'
      READ(*,'(A)') NAME
      K=INDEX(NAME,' ')
      IF (K.EQ.1) THEN
         PRINT *, 'Names may not start with a blank'
         GOTO 10
      ENDIF
      FILNAM = NAME(1:K-1)//'.rot' 
      OPEN(36,FILE=FILNAM,FORM='UNFORMATTED',STATUS = 'UNKNOWN')
      WRITE(36) 'G92RWF'

*      CALL CHKPLT

*   Perform machine- and installation-dependent setup
*
      CALL SETMC
*
*   Set up the physical constants
*
      CALL SETCON
*
*   Open, check, load data from, and close, the  .csl  file
*
      CALL SETCSLA(NAME,ncore_not_used)
*
*   Read the radial wave functions
*
      CALL GETHFD(NAME)

      WRITE(*,*)
      WRITE(*,*) ' Overlap matrix before rotation'
      WRITE(*,*)
      DO K = 1,NW
         DO I = 1,K
           IF (NAK(K).EQ.NAK(I)) THEN
             WRITE(*,1000) NP(K),NH(K),NP(I),NH(I),RINT(K,I,0)
           ENDIF
         ENDDO
      ENDDO
      WRITE(*,*)
      WRITE(*,*) ' Radius of orbitals '
      WRITE(*,*)
      DO K = 1,NW
         WRITE(*,1005) NP(K),NH(K),RINT(K,K,1),RINT(K,K,0)
      END DO

      write(*,*) 
      WRITE(*,*) '  s, p-, p , d-, d , f-, f, g-, g'
      WRITE(*,*) ' -1, 1 ,-2 , 2 ,-3 , 3 ,-4, 4 ,-5'
      WRITE(*,*)
      WRITE(*,*) ' Enter kappa for the orbital to be rotated'
      READ(*,*) NAKROT

      WRITE(*,*) ' How many do you want to skip ? '
      READ(*,*) ISKIP

      IC = 0
      DO J = 1,NW
         IF (NAK(J).EQ.NAKROT) THEN
            IFIRST=J
            IC = IC + 1
            IF (IC.EQ.ISKIP+1) GOTO 126
         ENDIF
      ENDDO

126   WRITE(*,1002) IFIRST
*
*   Rotate the orbital
*
      NFO = 1
      DO J = IFIRST+1,NW
         IF (NAK(J).EQ.NAKROT) THEN
            NFO = NFO + 1
            MFMAX = MAX(MF(IFIRST),MF(J))
            DO K = 1,MFMAX
               PF(K,IFIRST) = PF(K,IFIRST) + PF(K,J)
               QF(K,IFIRST) = QF(K,IFIRST) + QF(K,J)
            ENDDO
            PZ(IFIRST) = PZ(IFIRST) +  PZ(J) 
            MF(IFIRST) = MFMAX
         ENDIF
      ENDDO
      
      IF (NFO.EQ.1) THEN
        WRITE(*,*) ' Only one orbital of this symmetry... stop...'
        STOP
      ENDIf
      WRITE(*,1004) NP(IFIRST),NH(IFIRST),NFO-1
      WRITE(*,*)
      WRITE(*,*) ' Rotated  orbitals : '
      WRITE(*,*) ' ------------------- '
      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,*) ' Radius of rotated orbitals '
      WRITE(*,*)
*
*   Normalize the rotated orbital and save it
*
      QNORM = DSQRT(RINT(IFIRST,IFIRST,0))
      DO J = 1,MF(IFIRST)
         PF(J,IFIRST) = PF(J,IFIRST)/QNORM
         QF(J,IFIRST) = QF(J,IFIRST)/QNORM
      ENDDO
      PZ(IFIRST) = PZ(IFIRST)/QNORM 
*
      I = 1
   12 CONTINUE
C
C  *****  ORTHOGONALIZE TO INNER FUNCTIONS
C
      IM = I - 1
      DO II =1,IM
        IF (NAK(II) .NE. NAK(I)) GOTO 6
        PN = RINT(I,II,0)
        PNN = DSQRT(1.D0-PN*PN)
        MFMAX = MAX(MF(I),MF(II))
        DO J = 1,MFMAX
          PF(J,I) = (PF(J,I) - PN*PF(J,II))/PNN
          QF(J,I) = (QF(J,I) - PN*QF(J,II))/PNN
        ENDDO
        PZ(I) = (PZ(I) - PN*PZ(II))/PNN
        MF(I) = MFMAX
6     ENDDO

      PNN = RINT(I,I,0)   
      PNN = DSQRT(PNN)
      DO J = 1,MF(I)
        PF(J,I) = PF(J,I)/PNN
        QF(J,I) = QF(J,I)/PNN
      ENDDO
      PZ(I) = PZ(I)/PNN
      RAD = RINT(I,I,1)
      PNN = RINT(I,I,0)
      WRITE(*,1005) NP(I),NH(I),RAD,PNN

      WRITE(36) NP(I),NAK(I),E(I),MF(I)
      WRITE(36) PZ(I),(PF(J,I),J = 1,MF(I)),(QF(J,I),J = 1,MF(I))
      WRITE(36) (R(J),J = 1,MF(I))
      I = I + 1
9     IF ( I .LT. NW + 1) THEN
        GO TO 12
      ENDIF

      WRITE(*,*)
      WRITE(*,*) ' Overlap matrix after rotation'
      WRITE(*,*)
      DO K = 1,NW
         DO I = 1,K
           IF (NAK(K).EQ.NAK(I)) THEN
             WRITE(*,1000) NP(K),NH(K),NP(I),NH(I),RINT(K,I,0)
           ENDIF
         ENDDO
      ENDDO

1000  FORMAT(1H ,'<',I2,A2,'|',I2,A2,'> = ',F18.10)
1002  FORMAT(//,1H ,'The first orbital of this symmetry has # ',I3)
1004  FORMAT(//1H ,'The orbital ',I2,A2,' is  mixed (identical weights)
     : with the ',I4,//,' orbitals of the same symmetry. The latter
     : will be orthogonalized to the inner functions.',//)
1005  FORMAT(1H ,' Orbital  ',I2,A2,' with radius',F18.10,
     :  ' <|> = ',F18.10)

      END
************************************************************************
*                                                                      *
      SUBROUTINE GETHFD(NAME)
*                                                                      *
*   Interactively determines the data governing the HFS problem.       *
*                                                                      *
*   Call(s) to: [LIB92]: NUCPOT, RADGRD, SETQIC.                       *
*               [RCI92]: SETISO, SETRWF.                               *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 15 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = 600)
CGG      PARAMETER (NNNW = 120)
      INTEGER PNTRIQ
      LOGICAL GETYN,LFORDR,LTRANS,LVP,LSE,LNMS,LSMS,YES
      CHARACTER*24 NAME
*
      POINTER (PNTRPF,PF(NNNP,1)),(PNTRQF,QF(NNNP,1))
*
      COMMON/DECIDE/LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
     :      /DEF1/EMN,IONCTY,NELEC,Z
     :      /DEF2/C
     :      /DEF4/ACCY,NSCF,NSIC,NSOLV
     :      /DEF9/CVAC,PI
     :      /DEFAULT/NDEF
     :      /FOPARM/ICCUT
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /NPAR/PARM(2),NPARM
     :      /NPOT/ZZ(NNNP),NNUC
     :      /ORB1/E(NNNW),GAMA(NNNW)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
     :      /WFAC/WFACT
*
*   Open, check, load data from, and close the  .iso  file
*
      CALL SETISO ('isodata')
*
*   Determine the physical effects specifications
*
      IF (NDEF.NE.0) THEN
         PRINT *, 'The physical speed of light in'
         PRINT *, ' atomic units is',CVAC,';'
         PRINT *, ' revise this value?'
         YES = GETYN ()
         IF (YES) THEN
            PRINT *, 'Enter the revised value:'
            READ *,C
         ELSE
            C = CVAC
         ENDIF
      ELSE
         C = CVAC
      ENDIF
*
*   Determine the parameters controlling the radial grid
*
      IF (NPARM .EQ. 0) THEN
         RNT = EXP (-65.0D 00/16.0D 00) / Z
         H = 0.5D 00**4
         N = MIN (220,NNNP)
      ELSE
         RNT = 2.0D-06
         H = 5.0D-02
         N = NNNP
      ENDIF
      HP = 0.0D 00
      IF (NDEF.NE.0) THEN
         PRINT *, 'The default radial grid parameters'
         PRINT *, ' for this case are:'
         PRINT *, ' RNT = ',RNT,';'
         PRINT *, ' H = ',H,';'
         PRINT *, ' HP = ',HP,';'
         PRINT *, ' N = ',N,';'
         PRINT *, ' revise these values?'
         YES = GETYN ()
         IF (YES) THEN
            PRINT *, 'Enter RNT:'
            READ *, RNT
            PRINT *, 'Enter H:'
            READ *, H
            PRINT *, 'Enter HP:'
            READ *, HP
            PRINT *, 'Enter N:'
            READ *, N
         ENDIF
      ENDIF
*
*   ACCY is an estimate of the accuracy of the numerical procedures
*
      ACCY = H**6
*
*   Set up the coefficients for the numerical procedures
*
      CALL SETQIC
*
*   Generate the radial grid and all associated arrays
*
      CALL RADGRD
*
*   Load the radial wavefunctions
*
      CALL SETRWFA(TRIM(NAME)//'.w')
*
      RETURN
      END
