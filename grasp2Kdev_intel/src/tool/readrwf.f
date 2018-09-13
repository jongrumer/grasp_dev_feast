      PROGRAM readrwf
      IMPLICIT REAL*8 (A-H, O-Z)

!   Conver graspVU radial functions

      include 'parameters.def'
CGG      parameter (nnnp = 590)
CGG      parameter (nnn1 = nnnp+10)
CGG      PARAMETER (NNNW = 120)

      CHARACTER*6 G92RWF, InFile*80, OutFile*80
      DIMENSION PA(NNNP), QA(NNNP), RA(NNNP)

      PRINT *, '*******************************************************'
      PRINT *, '*'
      PRINT *, '* Convert graspVU radial function between two formats:'
      PRINT *, '* the UNFORMATTED and the FORMATTED (FORTRAN)'
      PRINT *, '*'
      PRINT *, '*******************************************************'
      PRINT *
   10 PRINT *, 'Input mode pls. (1: Unformatted to Formatted; 2: other)'
      READ (*,*, IOSTAT=ios) mode
      IF (ios .NE. 0 .OR. mode .GT. 2 .OR. mode .LT. 1) THEN
         PRINT *, 'Sorry, mode value not acceptable. Re-Enter pls'
         GOTO 10
      ENDIF
      PRINT *, 'Enter the input file pls'
      READ (*,'(A)') InFile
      PRINT *, 'Enter the output file pls:'
      READ (*,'(A)') OutFile
      PRINT *, 'Ok, I''ll do the rest'

      IF (mode .EQ. 1) THEN
         OPEN (23, FILE=InFile, FORM='UNFORMATTED',STATUS='OLD')
         OPEN (24, FILE=OutFile, STATUS='UNKNOWN')

         READ (23,IOSTAT = IOS) G92RWF
         IF ((IOS .NE. 0) .OR. (G92RWF .NE. 'G92RWF')) THEN
            WRITE (*,*) 'This is not a Radial WaveFunction File;'
            CLOSE (23)
            CLOSE (24)
            STOP
         ENDIF
         WRITE (24, '(A6)') G92RWF
         WRITE (6, '(A6)') G92RWF

         !.. Read orbitals until the end of the file

  100    CONTINUE

         READ (23, END = 200) NPY,NAKY,EY,MY
         WRITE (6, '( ''npy, naky, ey, my = '', 2I6, E25.16, I6 )') 
     &        NPY,NAKY,EY,MY
         IF (my .GT. NNNP) THEN
            WRITE (*,*) ' my = ', my, ' > ', 'NNNP = ', NNNP
            CLOSE (23)
            CLOSE (24)
            STOP 'Increase parameter NNNP and recompile'
         ENDIF
         READ (23) PZ,(PA(I),I = 1,MY),(QA(I),I = 1,MY)
         READ (23) (RA(I),I = 1,MY)

         WRITE (24, 1001) NPY,NAKY,EY,MY
         WRITE (24, 1002) PZ,(PA(I),I = 1,MY),(QA(I),I = 1,MY)
         WRITE (24, 1003) (RA(I),I = 1,MY)

         GOTO 100

  200    CONTINUE

      ELSE IF (mode .EQ. 2) THEN
         OPEN (23, FILE=InFile, STATUS='OLD')
         OPEN (24, FILE=OutFile, FORM='UNFORMATTED', STATUS='UNKNOWN')

         print *, 'ok 1'
         READ (23, '(A6)', IOSTAT = IOS) G92RWF
         print *, 'ok 2'
         IF ((IOS .NE. 0) .OR. (G92RWF .NE. 'G92RWF')) THEN
            WRITE (*,*) 'This is not a Radial WaveFunction File;'
            CLOSE (23)
            CLOSE (24)
            STOP
         ENDIF
         WRITE (24) G92RWF

         !.. Read orbitals until the end of the file

  300    CONTINUE

         print *, 'ok 3'
         READ (23, 1001, END = 400) NPY,NAKY,EY,MY
         IF (my .GT. NNNP) THEN
            WRITE (*,*) ' my = ', my, ' > ', 'NNNP = ', NNNP
            CLOSE (23)
            CLOSE (24)
            STOP 'Increase parameter NNNP and recompile'
         ENDIF
         print *, 'ok 4'
         READ (23,1002) PZ,(PA(I),I = 1,MY),(QA(I),I = 1,MY)
         print *, 'ok 5'
         READ (23,1003) (RA(I),I = 1,MY)
         print *, 'ok 6'

         WRITE (24) NPY,NAKY,EY,MY
         WRITE (24) PZ,(PA(I),I = 1,MY),(QA(I),I = 1,MY)
         WRITE (24) (RA(I),I = 1,MY)

         GOTO 300

  400    CONTINUE

      ENDIF

 1001 FORMAT (2I8, D25.16, I8)
 1002 FORMAT (D25.16)
 1003 FORMAT (D25.16)

      CLOSE (23)
      CLOSE (24)

      END
