*
************************************************************************
*                                                                      *
      PROGRAM ris_lsj
*                                                                      *
*     Program for displaying atomic istope shift  parameters           *
*     The program reads the LSJ classification                         *
*     file for labeling purposes.                                      *
*                                                                      *
*     Per Jonsson and Gediminas Gaigalas                               *    
*                                              November 2011           *
*                                                                      *
************************************************************************
      IMPLICIT NONE
      INTEGER, PARAMETER:: JMax = 22      ! max J value, see JFraction !
      INTEGER, PARAMETER:: ndim = 20000   ! max number of states
      INTEGER, PARAMETER:: maxFile = 1000 ! max number of files 
*
      CHARACTER(LEN=80) strInFile(maxFile), strFile
      CHARACTER*1 iaspa(ndim), PlusMinus(-1:1),ans    ! Parity
      CHARACTER*1 Lev_par(ndim),Lev_par_TOT(ndim),SWAPP                 ! Parity
      CHARACTER*4 iatjp(ndim), JFraction(1:2*Jmax+1)  ! J
      CHARACTER*4 Lev_J(ndim),Lev_J_TOT(ndim),SWAPJ                   ! J
      CHARACTER g92mix*6
      CHARACTER*64 string_CSF(ndim), string_PRN(ndim),SWAPCSF ! String in LSJ
      CHARACTER*64 string_CSF_TOT(ndim)
      CHARACTER*256 util_lbl_file, rms_file, output_file
      CHARACTER*32 DUMMY0
      CHARACTER*64 DUMMY
*
      INTEGER i, j, iargc, ios, ncountState, nFile, mFile
      INTEGER nelec, ncftot, nw, nvectot, nvecsiz, nblock, jblock
      INTEGER nb, ncfblk, nevblk, iiatjp, iiaspa
      INTEGER K, Iprint, IMaxCount, nsort, IMaxCount_tot
      INTEGER ivec(ndim), indx(ndim), Lev_POS(ndim), izero
      INTEGER MAX_STRING_LENGTH,MAX_STRING_LENGTH_TOT,N
*
      DOUBLE PRECISION eav, eval(ndim), evec, RLev_ENER(ndim),ZERO
      DOUBLE PRECISION K1(ndim),K2(ndim),K3(ndim),SWAP1,SWAP2,SWAP3,
     :                 SWAP4,N1(ndim),N2(ndim),N3(ndim),DENS(ndim),
     :                 RLev_ENER_TOT(ndim),K3_TOT(ndim),N3_TOT(ndim),
     :                 DENS_TOT(ndim)
*
      COMMON/JJ2LSJ/ Lev_POS,Lev_J,Lev_Par,RLev_ENER,string_CSF,
     :                IMaxCount , MAX_STRING_LENGTH
      COMMON/RMS/K1,K2,K3,DENS,N1,N2,N3
*
      DATA PlusMinus/'-', ' ', '+'/
      DATA JFraction/'  0 ', ' 1/2', '  1 ', ' 3/2', '  2 ', ' 5/2',
     &               '  3 ', ' 7/2', '  4 ', ' 9/2', '  5 ', '11/2',
     &               '  6 ', '13/2', '  7 ', '15/2', '  8 ', '17/2',
     &               '  9 ', '19/2', ' 10 ',
     &                       '21/2', ' 11 ', '23/2', ' 12 ', '25/2',
     &               ' 13 ', '27/2', ' 14 ', '29/2', ' 15 ', '31/2',
     &               ' 16 ', '33/2', ' 17 ', '35/2', ' 18 ', '37/2',
     &               ' 19 ', '39/2', ' 20 ', '41/2', ' 21 ', '43/2',
     &               ' 22 '/
*
*
*     Opens the file  *.lsj.lbl
*
      WRITE(*,*)
      WRITE(*,*) 'RIS_LSJ'
      WRITE(*,*) 'This program prints output from the ris program '
      WRITE(*,*) 'using LSJ lables. Output can be energy sorted'
      WRITE(*,*) 'Input files: name.(c)i, name.lsj.lbl'
      WRITE(*,*) 'Output file: name.(c)ilsj'
      WRITE(*,*)

      IMaxCount_tot = 0
      MAX_STRING_LENGTH_TOT = 0

      DO 
        WRITE(*,*) 'Name of the state'
        READ(*,*) strFile 
        K = INDEX(strFile,' ')
        WRITE(*,*) 'Isotope shift data from CI calculation?'
        READ(*,*) ans
        IF ((ans.eq.'y').or.(ans.eq.'Y')) THEN
           rms_file = strFile(1:K-1)//'.ci'
        ELSE
           rms_file = strFile(1:K-1)//'.i'
        END IF
        OPEN (30, FILE = rms_file, FORM = 'FORMATTED',
     &   STATUS = 'OLD', IOSTAT = IOS)
        IF (IOS .NE. 0) THEN
           WRITE (0,*) 'Failed to open file ', TRIM(rms_file)
           CLOSE (30)
           STOP
        END IF

        util_lbl_file = strFile(1:K-1)//'.lsj.lbl'
        OPEN (31, FILE = util_lbl_file, FORM = 'FORMATTED', 
     &     STATUS = 'OLD', IOSTAT = IOS)
        IF (IOS .NE. 0) THEN
           WRITE (0,*) 'Failed to open file ', TRIM(util_lbl_file)
           CLOSE (31)
           STOP
        END IF
        IMaxCount = 0
        CALL LDLBL
        CALL READRMS

        DO J = 1,IMaxCount
           RLev_ENER_TOT(J+IMaxCount_tot) = RLev_ENER(J)
           String_CSF_TOT(J+IMaxCount_tot) = String_CSF(J)
           Lev_J_TOT(J+IMaxCount_tot) = Lev_J(J)
           Lev_Par_TOT(J+IMaxCount_tot) = Lev_Par(J)
           K3_TOT(J+IMaxCount_tot) = K3(J)
           N3_TOT(J+IMaxCount_tot) = N3(J)
           DENS_TOT(J+IMaxCount_tot) = DENS(J)
        END DO
        IMaxCount_tot = IMaxCount_tot + IMaxCount
        IF (MAX_STRING_LENGTH.GT.MAX_STRING_LENGTH_TOT) THEN
           MAX_STRING_LENGTH_TOT = MAX_STRING_LENGTH
        END IF

        CLOSE (30)
        CLOSE (31)

CPJ        WRITE(*,*) 'More files?'
CPJ        READ(*,*) ans
CPJ        IF ((ans.eq.'n').or.(ans.eq.'N')) GOTO 73
           GOTO 73
      END DO
   73 CONTINUE

      IF ((ans.eq.'y').or.(ans.eq.'Y')) THEN
         OPEN (80, FILE = trim(strFile)//'.cilsj', FORM = 'FORMATTED',
     &     STATUS = 'UNKNOWN')
      ELSE
         OPEN (80, FILE = trim(strFile)//'.ilsj', FORM = 'FORMATTED',
     &     STATUS = 'UNKNOWN')
      END IF

      WRITE(*,*) 'Energy sorted output? '
      READ(*,*) ans
      IF ((ans.eq.'y').or.(ans.eq.'Y')) THEN
        nsort = 1
      ELSE
        nsort = 0
      END IF

      DUMMY0 = '                                '
      DUMMY = DUMMY0//DUMMY0
      N = MAX(1,MAX_STRING_LENGTH_TOT-4)
      WRITE(80,402) 'Energy',DUMMY(1:N),'State',
     &   '    J   P    K_NMS (a.u.)   K_SMS (a.u.)     Dens (a.u.)'

      IF (nsort.eq.1) THEN
        DO J = 2,IMaxCount_tot
          SWAP1 = RLev_ENER_TOT(J)
          SWAP2 = K3_TOT(J)
          SWAP3 = N3_TOT(J)
          SWAP4 = DENS_TOT(J)
          SWAPP = Lev_Par_TOT(J)
          SWAPJ = Lev_J_TOT(J)
          SWAPCSF = string_CSF_TOT(J)
          DO I = J-1,1,-1
            IF (RLev_ENER_TOT(I).LE.SWAP1) GOTO 10
            RLev_ENER_TOT(I+1) = RLev_ENER_TOT(I)
            K3_TOT(I+1) = K3_TOT(I)
            N3_TOT(I+1) = N3_TOT(I)
            DENS_TOT(I+1) = DENS_TOT(I)
            Lev_Par_TOT(I+1) = Lev_Par_TOT(I)
            Lev_J_TOT(I+1) = Lev_J_TOT(I)
            string_CSF_TOT(I+1) = string_CSF_TOT(I)
          END DO
          I = 0
   10     RLev_ENER_TOT(I+1) = SWAP1
          K3_TOT(I+1) = SWAP2
          N3_TOT(I+1) = SWAP3
          DENS_TOT(I+1) = SWAP4
          Lev_Par_TOT(I+1) = SWAPP
          Lev_J_TOT(I+1) = SWAPJ
          string_CSF_TOT(I+1) = SWAPCSF
        END DO
      END IF

      DO I = 1,IMaxCount_tot
         WRITE(80,405) RLev_ENER_TOT(I),
     &   string_CSF_TOT(I)(1:MAX_STRING_LENGTH_TOT),Lev_J_TOT(I),
     &   Lev_Par_TOT(I),N3_TOT(I),K3_TOT(I),DENS_TOT(I)
      END DO

      CLOSE (80)

      STOP 

  402 FORMAT (//9X,A,A,A,A)
C  402 FORMAT (//' Interaction constants:'//
C     :' Energy      J Parity ',1X,'A (MHz)',6X,'B (MHz)',8X,
C     : 'g_J            state'/)
  403 FORMAT (1X,1I3,5X,2A4,1P,2D13.3,D16.6,5X,A)
  404 FORMAT (1X,1I3,5X,2A4,1P,2D13.3,D16.6,5X,A,F14.7)
  405 FORMAT (1X,F14.7,2X,A,2A4,1P,3D16.7)


      END

*
************************************************************************
*                                                                      *
      SUBROUTINE READRMS
*                                                                      *
*     Open, check and load data from the ris3 file                      *
*                                                                      *
*     Calls:                                                           *
*                                                                      *
*     Written by P. J\"onsson,                                         *
*     Malmo                                                 Aug 2011   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      INTEGER, PARAMETER :: ndim = 20000          ! max number of states

      CHARACTER*1 Lev_par(ndim)
      CHARACTER*4 Lev_J(ndim)
      CHARACTER*64 string_CSF(ndim)
*
      INTEGER Lev_POS(ndim)

      CHARACTER*120 RECORD
      DOUBLE PRECISION K1(ndim),K2(ndim),K3(ndim),DENS(ndim)
      DOUBLE PRECISION N1(ndim),N2(ndim),N3(ndim)
      DOUBLE PRECISION RLev_ENER(ndim)
*
      COMMON/JJ2LSJ/ Lev_POS,Lev_J,Lev_Par,RLev_ENER,string_CSF,
     :                IMaxCount , MAX_STRING_LENGTH
      COMMON/RMS/K1,K2,K3,DENS,N1,N2,N3

* Position yourself at the correct place in the file

      DO
         READ (30,'(A)') RECORD
CPJ         WRITE(*,'(A)') RECORD
         IF (RECORD(19:24).EQ.'Normal') GOTO 10
      END DO
   10 CONTINUE

* Now begin to read the specific mass shift  data

      WRITE(*,*) IMaxCount
      DO I = 1,IMaxCount
         READ (30,'(A)') RECORD
         READ (30,'(A)') RECORD
         READ (30,'(20X,3D20.10)') N1(I),N2(I),N3(I)
!        WRITE(*,*) N1(I),N2(I),N3(I)
         READ (30,'(A)') RECORD
      END DO

      DO I = 1,3
         READ (30,'(A)') RECORD
!        write(*,*) RECORD
      END DO

      WRITE(*,*) IMaxCount
      DO I = 1,IMaxCount
!        write(*,*) 'I=',I
         READ (30,'(A)') RECORD
!        write(*,*) 'APA',RECORD
         READ (30,'(A)') RECORD
!        write(*,*) 'BANAN',RECORD
         READ (30,'(20X,3D20.10)') K1(I),K2(I),K3(I)
!           WRITE(*,*) K1(I),K2(I),K3(I)         
         READ (30,'(A)') RECORD
         write(*,*) RECORD
      END DO

      DO I = 1,6
         READ (30,'(A)') RECORD
      END DO

      DO I = 1,IMaxCount
         READ (30,'(20X,D20.10)') DENS(I)
CPJ         WRITE(*,*) DENS(I)
      END DO

!      DO I = 1,4
!         READ (30,'(A)') RECORD
!      END DO

!      DO I = 1,IMaxCount
!         READ (30,'(A)') RECORD
!         READ (30,'(A)') RECORD
!         READ (30,'(20X,3D20.10)') N1(I),N2(I),N3(I)
!         WRITE(*,*) N1(I),N2(I),N3(I)         
!         READ (30,'(A)') RECORD
!      END DO

      RETURN
      END

*
************************************************************************
*                                                                      *
      SUBROUTINE LDLBL
*                                                                      *
*     Open, check and load data from the  .lsj.lbl   file of the       *
*     inital state.                                                    *
*                                                                      *
*     Calls:                                                           *
*                                                                      *
*     Written by G. Gaigalas,                                          *
*     NIST                                                  May 2011   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      INTEGER, PARAMETER :: ndim = 20000          ! max number of states
*
      CHARACTER*1 Lev_par(ndim)
      CHARACTER*4 Lev_J(ndim)
      CHARACTER*15 RECORD
      CHARACTER*64 string_CSF(ndim)
*
      INTEGER IOS, ITEST, Lev_POS(ndim)
      INTEGER MAX_STRING_LENGTH
      REAL WEIGHTS
      DOUBLE PRECISION RLev_ENER(ndim)
*
      COMMON/JJ2LSJ/ Lev_POS,Lev_J,Lev_Par,RLev_ENER,string_CSF,
     :                IMaxCount, MAX_STRING_LENGTH
*
      MAX_STRING_LENGTH = 0 

      READ (31,'(1A15)',IOSTAT = IOS) RECORD
      ICount = 0
      IF (IOS .NE. 0) GO TO 1
      ICount = 1
      READ (31,'(1X,I3,1X,A4,5X,A1,8X,F16.9)',IOSTAT = IOS)
     :  Lev_Pos(ICount),Lev_J(ICount),Lev_Par(ICount),
     :  RLev_ENER(ICount)
      IF (IOS .NE. 0) GO TO 1
*
      READ (31,'(7X,F12.8,17X,A)') WEIGHTS,string_CSF(ICount)
      K = INDEX(string_CSF(ICount),' ')
      IF (K.GT.MAX_STRING_LENGTH) MAX_STRING_LENGTH = K
               
*
    2 READ (31,'(1X,I2)',IOSTAT = IOS) ITEST
      IF (IOS .NE. 0) GO TO 1
      IF (ITEST .EQ. 0) GO TO 2
      BACKSPACE 31
      ICount = ICount + 1
      READ (31,'(1X,I3,1X,A4,5X,A1,8X,F16.9)',IOSTAT = IOS)
     :  Lev_Pos(ICount),Lev_J(ICount),Lev_Par(ICount),
     :  RLev_ENER(ICount)
      IF (IOS .NE. 0) GO TO 1
      READ (31,'(7X,F12.8,17X,A)') WEIGHTS,string_CSF(ICount)
      K = INDEX(string_CSF(ICount),' ')
      IF (K.GT.MAX_STRING_LENGTH) MAX_STRING_LENGTH = K
      GO TO 2
    1 CONTINUE
      IMaxCount = ICount + IMaxCount
      WRITE(*,*) IMaxCount,1
      RETURN
      END
