*
************************************************************************
*                                                                      *
      PROGRAM rlevelsname
*                                                                      *
*     Purpose:                                                         *
*     Extract energy levels from rscf/rci output files (.m, .cm file)  *
*     and from jj2lsj output file (*.lsj.lbl).                         *
*     And then:                                                        *
*         (1) sort                                                     *
*         (2) difference with the lowest                               *
*         (3) difference with the nearby lower                         *
*     Files names are provided by user.                                *
*     Usage:                                                           *
*           $ grlevels file1 file2 file3 ...                           *
*           or                                                         *
*           $ grlevels                                                 *
*           file1                                                      *
*           file2                                                      *
*           file3                                                      *
*           ... (return to terminate)                                  *
*                                                                      *
*     Calls: LDLBL, indexS.                                            *
*                                                                      *
*     Xinghong He  98-10-16                                            *
*                                                                      *
*     Rewritten by  G. Gaigalas                                        *    
*     for LSJ calssification of levels                                 *    
*     NIST                                                 May 2011    *
*     Adapted by Per Jonsson for producing latex tables November 2014  *      
*                                                                      *
************************************************************************
      IMPLICIT NONE
      include 'parameters.def'
      INTEGER, PARAMETER:: JMax = 22      ! max J value, see JFraction !
      INTEGER, PARAMETER:: ndim = 20000   ! max number of states
      INTEGER, PARAMETER:: maxFile = 1000 ! max number of files 
      DOUBLE PRECISION, PARAMETER:: Rydberg = 109737.31534D0
*
      CHARACTER(LEN=20) strInFile(maxFile), strFile, strFile2
      CHARACTER*1 iaspa(ndim), PlusMinus(-1:1)        ! Parity
      CHARACTER*1 Lev_par(ndim)                       ! Parity
      CHARACTER*4 iatjp(ndim), JFraction(1:2*Jmax+1)  ! J
      CHARACTER*4 Lev_J(ndim)                         ! J
      CHARACTER g92mix*6
      CHARACTER*64 string_CSF(ndim), string_PRN(ndim), latexstring ! String in LSJ
      CHARACTER*20 string_file(ndim)
      CHARACTER util_lbl_file*256
*
      INTEGER i, j, iargc, ios, ncountState, nFile, mFile, nfound
      INTEGER nelec, ncftot, nw, nvectot, nvecsiz, nblock, jblock
      INTEGER nb, ncfblk, nevblk, iiatjp, iiaspa
      INTEGER K, Iprint, IMaxCount
      INTEGER ivec(ndim), indx(ndim), Lev_POS(ndim)
*
      DOUBLE PRECISION eav, eval(ndim), evec, RLev_ENER(ndim),ZERO
*
      COMMON/JJ2LSJ/ Lev_POS,Lev_J,Lev_Par,RLev_ENER,string_CSF,
     :                IMaxCount 
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
      OPEN(UNIT=9,file='energylabel',status='unknown',form='formatted')

      mFile = iargc()
      string_PRN = ''
      Iprint = 0
      ZERO = 0.0D00
      IF (mFile .EQ. 0) THEN                 ! Get file names interactively
         WRITE (0,*)'  You can also use command-line option:'
         WRITE (0,*)'    %grlevels file1 file2 (wild cards allowed)...'
         WRITE (0,*)'  Now, carry on'
         WRITE (0,*)
         WRITE (0,*)'Type the input file name, one for each line',
     &            ' (NULL to terminate)'
         WRITE (0,*)
         i = 0
         DO  ! Don't know the exact number of files.
            WRITE (0,'(A12)', ADVANCE='NO') 'File name ? '
            READ (5, '(A20)') strFile
            strFile = ADJUSTL (strFile)
            IF (LEN_TRIM (strFile) .GT. 0) THEN ! a valid input
               i = i + 1
               IF (i .GT. maxFile) THEN      !  impose an upper limit
                  WRITE (0,*) 'Too many files opened. processing first '
     &                       , i-1, ' files.'
                  EXIT
               ENDIF
               strInFile(i) = strFile
            ELSE
               EXIT
            ENDIF
         ENDDO
         mFile = i 
      ELSEIF (mFile .GT. 0 .AND. mFile .LE. maxFile) THEN
         DO i = 1, mFile
            CALL getarg (i, strInFile(i))
         ENDDO
      ELSE
         WRITE (0,*) 'More than ', maxFile, ' files entered,',
     &            ' modify parameter maxFile'
      ENDIF
*
*     Open files, read energies, concatenate to a single place
*     Open mix file, check header
*
      ncountState = 0
      IMaxCount = 0
      DO nFile = 1, mFile
         strFile = strInFile(nFile)
         OPEN (3, FILE = strFile, FORM = 'UNFORMATTED', STATUS = 'OLD'
     &      , IOSTAT = IOS)
         IF (IOS .NE. 0) THEN
           WRITE (0,*) 'Failed to open file "', 
     &          strFile(1:LEN_TRIM (strFile)), '", skipping...'
            CLOSE (3)
            CYCLE
         ENDIF
         READ (3) g92mix
         IF (g92mix .NE. 'G92MIX') THEN
            WRITE (0,*) 'Not a mixing coefficient file, skipping "',
     &                strFile(1:LEN_TRIM (strFile)), '"'
            CLOSE (3)
            CYCLE
         ENDIF

* Now trim the file string

         strFile2 = strFile
         nfound = 0
         DO i = 20,1,-1
            if (strFile2(i:i).eq.'/') then
               nfound = i
               exit
            end if 
         END DO
         if (nfound.ne.0) then
            do i = 1,nfound
               strFile2(i:i) = ' '
            end do
         end if
         nfound = 0
         DO i = 20,3,-1
            if ((strFile2(i-2:i).eq.'.cm').OR.
     &         (strFile2(i-2:i).eq.'.m ')) 
     &         then
               nfound = i
               exit
            end if 
         END DO
         if (nfound.ne.0) then
            do i = 1,3
               strFile2(i-3+nfound:i-3+nfound) = ' '
            end do
         end if

         strFile2 = adjustl(strFile2)

         READ (3) nelec, ncftot, nw, nvectot, nvecsiz, nblock
*
         WRITE(9,*) 'nblock = ', nblock, '  ncftot = ', ncftot, 
     &                '  nw = ', nw, '  nelec = ', nelec

         DO jblock = 1, nblock
            READ (3) nb, ncfblk, nevblk, iiatjp, iiaspa
            IF (jblock .NE. nb) THEN
*
*     This error can occur anywhere and therefore cannot 
*     be simply skipped - stop instead.
*
               WRITE (0,*) 'jblock .NE. nb, stopping...'
               CLOSE (3)
               STOP
            ENDIF
            IF (nevblk .LE. 0) CYCLE
            READ (3) (ivec(i+ncountState), i = 1, nevblk)
            READ (3) eav, (eval(i+ncountState), i = 1, nevblk)
            READ (3) (evec, i = 1, ncfblk*nevblk)
*
*     Assign J and parity to every individual state
*     Also add the average energy (back) to energy
*
            DO i = 1, nevblk
               iatjp(i+ncountState) = JFraction(iiatjp)
               iaspa(i+ncountState) = PlusMinus(iiaspa)
               eval(i+ncountState) = eval(i+ncountState) + eav
               string_file(i+ncountState) = strFile2
            END DO
*
*     Update ncountState
*
            ncountState = ncountState + nevblk
         END DO
         CLOSE (3)
*
*     Opens the file  *.lsj.lbl
*
         K = INDEX(strFile,' ')
         if(strFile(K-2:K-1) .EQ. '.m') then
            util_lbl_file = strFile(1:K-3)//'.lsj.lbl'
         else if(strFile(K-3:K-1) .EQ. '.cm') then
            util_lbl_file = strFile(1:K-4)//'.lsj.lbl'
         end if
         OPEN (31, FILE = util_lbl_file, FORM = 'FORMATTED', 
     &      STATUS = 'OLD', IOSTAT = IOS)
         IF (IOS .NE. 0) THEN
CGG            WRITE (0,*) 'Failed to open file "', 
CGG     &      util_lbl_file(1:LEN_TRIM (util_lbl_file)), '", skipping...'
            CLOSE (31)
            CYCLE
         END IF
         Iprint = 1
*
*     Defines the LSJ string for the levels
*
         CALL LDLBL
         CLOSE (31)
         do i = 1,ncountState
            do j = 1,IMaxCount
               if(DABS(DABS(eval(i))-DABS(RLev_ENER(j))) 
     :                           .LT. 0.00000001) then
                   string_PRN(i) =  string_CSF(j)
                   exit
               end if
            end do
         end do
      END DO
      CALL indexS (ncountState, eval, .FALSE., indx)
*
*     The output of the levels
*
      WRITE (9,*)
      WRITE (9,1) 
      WRITE (9,2) Rydberg
      if(Iprint .eq. 1) then
         WRITE (9,5)
         WRITE (9,*) 'No - Serial number of the state; ',
     &               'Pos - Position of the state within the '
         WRITE (9,*) 'J/P block;' 
         WRITE (9,5)
         WRITE (9,*) 'No Pos  J Parity Energy Total    Levels',
     &               '   File          Configuration'
         WRITE (9,*) '                     (a.u.)      (cm^-1)'
         WRITE (9,5)
         j = 1
         i = indx(j)
         call latexconvert(string_PRN(i),latexstring)
         latexstring = '$'//trim(latexstring)//'_{'//iatjp(i)
         if (iaspa(i).eq.'-') then
            latexstring = trim(latexstring)//'}^o$'
         else
            latexstring = trim(latexstring)//'}$'
         end if

         WRITE (9,3) j,ivec(i),iatjp(i),iaspa(i),eval(i),
     &      ZERO,
     &      string_file(i),Trim(latexstring)
!     &      string_file(i),Trim(string_PRN(i))
         DO j = 2, ncountState
            i = indx(j)
            call latexconvert(string_PRN(i),latexstring)
            latexstring = '$'//trim(latexstring)//'_{'//iatjp(i)
            if (iaspa(i).eq.'-') then
               latexstring = trim(latexstring)//'}^o$'
            else
               latexstring = trim(latexstring)//'}$'
            end if
            WRITE (9,3) j,ivec(i),iatjp(i),iaspa(i),eval(i),
     &      (eval(i)-eval(indx(1)))*Rydberg*2,
     &      string_file(i),Trim(latexstring) 
!     &      string_file(i),Trim(string_PRN(i)) 
         END DO
!         WRITE (6,5)
      ELSE
         WRITE (9,*) 'No - Serial number of the state; ',
     &               'Pos - Position of the state within the '
         WRITE (9,*) 'J/P block;'
         WRITE (9,5)
         WRITE (9,*) 'No Pos  J Parity Energy Total    Levels',
     &               '   File          Configuration'
         WRITE (9,*) '                     (a.u.)      (cm^-1)'
         WRITE (9,5)
         j = 1
         i = indx(j)
         write(latexstring,'(a1,i3,a3,a4,3a1)')  '$',j,'~J=',
     &            iatjp(i),'^',iaspa(i),'$'
         WRITE (9,3)
         WRITE (9,3) j,ivec(i),iatjp(i),iaspa(i),eval(i),ZERO,
     &              string_file(i),trim(latexstring)
         DO j = 2, ncountState
            i = indx(j)
            write(latexstring,'(a1,i3,a3,a4,3a1)')  '$',j,'~J=',
     &           iatjp(i),'^',iaspa(i),'$'
            WRITE (9,3) j,ivec(i),iatjp(i),iaspa(i),eval(i),
     &      (eval(i)-eval(indx(1)))*Rydberg*2,
     &      string_file(i),trim(latexstring)
         END DO
!         WRITE (6,6)
      END IF

      write(*,*) ' Energy label data written to file'
      write(*,*) ' energylabel'

   1  FORMAT (' Energy levels for ...')
   2  FORMAT (' Rydberg constant is ',F14.5)
   3  FORMAT (2I3,1X,A4,1x,A1,2X,F14.7,F12.2,2X,A,A)
   5  FORMAT ('---------------------------------------------',
     &        '---------------------------------------------')
   6  FORMAT ('------------------------------------------',
     &        '-------------------------------')
      STOP
      END
*
************************************************************************
*                                                                      *
      SUBROUTINE indexS(n,a,ldown,indx)
*                                                                      *
*     Sort out the order of array a and store the index in indx        *
*                                                        (a pointer)   *
*     The input array a is unchanged  written in the bases of UpDown   *
*                                                                      *
*     !$Id: rlevels.f,v 1.2 2003/10/02 07:56:22 per Exp $              *
*     $Log: rlevels.f,v $                                              *
*     Revision 1.2  2003/10/02 07:56:22  per                           *
*     *** empty log message ***                                        *
*                                                                      *
*     Revision 1.1.1.1  2003/01/04 21:45:39  georgio                   *
*                                                                      *
************************************************************************
      IMPLICIT NONE
      CHARACTER*(*)    RCSID
      PARAMETER        ( RCSID
     & ='$Id: rlevels.f,v 1.2 2003/10/02 07:56:22 per Exp $'
     & )
      LOGICAL          ldown              ! .TRUE. then Big ---> Small
      INTEGER          n, indx(n)
      DOUBLE PRECISION a(n)
      INTEGER          i, j, ipos, jpos, jhere
      DOUBLE PRECISION aimx
!
!     Initialize the index array
      DO i = 1, n
         indx(i) = i
      ENDDO
      IF (ldown) THEN
         DO i = 1, n
            ipos = indx(i)
            aimx = a(ipos)
            jhere = i
            DO j = i+1, n
               jpos = indx(j)
               IF(a(jpos) .GT. aimx) THEN
                  aimx = a(jpos)
                  jhere = j
               ENDIF
            ENDDO
            indx(i) = indx(jhere)
            indx(jhere) = ipos
         ENDDO
      ELSE
         DO i = 1, n
            ipos = indx(i)
            aimx = a(ipos)
            jhere = i
            DO j = i+1, n
               jpos = indx(j)
               IF(a(jpos) .LT. aimx) THEN
                  aimx = a(jpos)
                  jhere = j
               ENDIF
            ENDDO
            indx(i) = indx(jhere)
            indx(jhere) = ipos
         ENDDO
      ENDIF
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
      REAL WEIGHTS
      DOUBLE PRECISION RLev_ENER(ndim)
*
      COMMON/JJ2LSJ/ Lev_POS,Lev_J,Lev_Par,RLev_ENER,string_CSF,
     :                IMaxCount 
*
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
*
    2 READ (31,'(1X,I2)',IOSTAT = IOS),ITEST
      IF (IOS .NE. 0) GO TO 1
      IF (ITEST .EQ. 0) GO TO 2
      BACKSPACE 31
      ICount = ICount + 1
      READ (31,'(1X,I3,1X,A4,5X,A1,8X,F16.9)',IOSTAT = IOS)
     :  Lev_Pos(ICount),Lev_J(ICount),Lev_Par(ICount),
     :  RLev_ENER(ICount)
      IF (IOS .NE. 0) GO TO 1
      READ (31,'(7X,F12.8,17X,A)') WEIGHTS,string_CSF(ICount)
      GO TO 2
    1 CONTINUE
      IMaxCount = ICount + IMaxCount
      RETURN
      END

******************************************************************

      subroutine latexconvert(labelstring,latexstring)

! This subroutine converts a label string to latex
! It is basically the same routine as in renergytable.f90      
! Per Jonssson, Malmo University, November 2014

      implicit none
      integer :: i,j,k,l,ncase
      character(len=64) :: labelstring, dummystring, latexstring
      character(len=1) :: char1, char2, char3

      do i = 1,61

!  Replace (n) with ^n

         if ((labelstring(i:i).eq.'(').and.
     &      (labelstring(i+2:i+2).eq.')')) then
            labelstring(i:i) = '^'
            labelstring(i+2:i+2) = ' '
         end if
      end do

      do i = 1,61

!  Replace . with \,

         if (labelstring(i:i).eq.'.') then
            dummystring = labelstring
            labelstring(1:i-1) = dummystring(1:i-1)
            labelstring(i:i) = '\'
            labelstring(i+1:i+1) = ','
            labelstring(i+2:64) = dummystring(i+1:62) 
         end if
      end do

      do i = 1,61

!  Replace _ with ~

         if (labelstring(i:i).eq.'_') labelstring(i:i) = '~'
      end do

!  If integer1 and S, P, D, F, G, H, I, K, L, M, N and integer2 replace with (^integer1_integer2S), (^integer1_integer2P), etc

      do l = 1,15
         ncase = 0
         do i = 1,61
            do j = 48,57
               do k = 48,57 
                  char1 = labelstring(i:i)
                  char2 = labelstring(i+1:i+1)
                  char3 = labelstring(i+2:i+2)
                  if ((ichar(char1).eq.j).and.(ichar(char3).eq.k).and.
     &               ((char2.ne.'~').and.(char2.ne.' ').and.
     &               (char2.ne.'_'))) then
                     dummystring = labelstring
                     labelstring(1:i-1) = dummystring(1:i-1)
                     labelstring(i:i+6) = 
     &                    '(^'//char1//'_'//char3//char2//')'
                     labelstring(i+7:64) = dummystring(i+3:60)
                     ncase = ncase + 1
                  end if
               end do
            end do
            if (ncase.eq.1) exit
         end do

      end do

!  If integer1 and S, P, D, F, G, H, I, K, L, M, N and not integer2 replace with ^integer1S, ^integer1P, etc 

      do i = 1,61
!  
         if (labelstring(i:i).eq.'~') then
            dummystring = labelstring
            labelstring(1:i) = dummystring(1:i)
            labelstring(i+1:i+1) = '^'
            labelstring(i+2:64) = dummystring(i+1:62)
         end if
      end do

      latexstring = trim(labelstring)

      return
      end      


