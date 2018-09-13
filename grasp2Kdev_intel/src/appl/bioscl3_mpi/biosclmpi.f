************************************************************************
*                                                                      *
      PROGRAM BIOSCL
*                                                                      *
*   This program calculates the transition parameters for a            *
*   transition between an initial and a final state                    *
*   The program assumes that the radial orbitals of the initial        *
*   and final state have been transformed by the BIOTRA program        *
*   as to become biorthonormal, in which case the normal Racah         *
*   algebra can be used.                                               *
*                                                                      *
*   Written by Per Jonsson,   Department of Computer Science           *
*                             Vanderbilt University, USA               *
*                                                                      *
*   MPI -- Rasa Matulioniene, July 2000                                *
*   Modified by Gediminas Gaigalas for new spin-angular integration.   *
*   Modify by G. Gaigalas                                12 Dec 2014   *
************************************************************************
*
      LOGICAL GETYN,YES,LDBPR
      CHARACTER*24 NAME(2)
      CHARACTER*128 FULLNAME(2)  ! includes pathname
      DOUBLE PRECISION CUTOFF

      COMMON/DEFAULT/NDEF,NDUMP,INPCI
      COMMON/DEBUGR/LDBPR(30), CUTOFF
CGG  lbl
      COMMON/JJ2LSJ3/ tmpdir,startdir,idstring


      INCLUDE 'mpif.h'
      INTEGER myid, nprocs, ierr, lenhost
      COMMON /mpi/ myid, nprocs, ierr
      CHARACTER host*(MPI_MAX_PROCESSOR_NAME), idstring*3

!Things for timing
      INTEGER ncount1
 
! cpath uses
      INTEGER       lenperm, lentmp
      CHARACTER*128 startdir, permdir, tmpdir
!-----------------------------------------------------------------------
 
*=======================================================================
*  Start mpi --- get processor info: myid, nprocs, host name; and print
*=======================================================================        
      CALL startmpi2 (myid, nprocs, host, lenhost, ncount1,
     &                startdir, permdir, tmpdir, 'RTRANSITION_MPI')
      WRITE (idstring, '(I3.3)') myid
      print*, tmpdir , ' = tmpdir'
 
      lenperm = LEN_TRIM (permdir)
      lentmp = LEN_TRIM (tmpdir)                                                

*=======================================================================        

      NTEST = 1001

      IF (myid .EQ. 0) THEN
          PRINT *
          PRINT *, ' Default settings?'
          YES = GETYN ()
          YES = .TRUE.
          IF (YES) THEN
            NDEF = 0
            NDUMP = 1
          ELSE
            NDEF = 1
            PRINT *, ' Dump angular data to file?'
            YES = GETYN ()
            IF (YES) THEN
              NDUMP = 1
            ELSE
              NDUMP = 0
            ENDIF
          ENDIF
          PRINT *, ' Input from a CI calculation?'
          YES = GETYN ()
          PRINT *
          IF (YES) THEN
            INPCI = 0
          ELSE
            INPCI = 1
          ENDIF
      ENDIF   !myid=0
      CALL MPI_Bcast (NDEF,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)   
      CALL MPI_Bcast (NDUMP,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)   
      CALL MPI_Bcast (INPCI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)   
*Rasa -- start
*     LDBPR = .F.
*     PRINT *, ' Generate debug output?'
*     YES = GETYN ()
*     PRINT *
*     IF (YES) THEN
*       LDBPR(18) = .T.
*       PRINT *, ' Enter the cutoff'
*       read *, cutoff
*     ENDIF
*Rasa -- end
*
*   Perform machine- and installation-dependent setup
*
      CALL SETMC
*
*   Set up the physical constants
*
      CALL SETCON
*
*   Obtain the names of the initial and final state files
*
      IF (myid .EQ. 0) CALL FNAME(NAME)
      CALL MPI_Bcast (NAME,48,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)   
*R    print *, 'DBG: myid=',myid, ' name(1)=', name(1)
*R    print *, 'DBG: myid=',myid, ' name(2)=', name(2)
      FULLNAME(1) = TRIM(startdir) // '/' // TRIM(NAME(1))
      FULLNAME(2) = TRIM(startdir) // '/' // TRIM(NAME(2))
*R    print *, 'DBG: myid=',myid, ' startdir =', startdir
*R    print *, 'DBG: myid=',myid, ' fullname(1)=', fullname(1)
*R    print *, 'DBG: myid=',myid, ' fullname(2)=', fullname(2)
*
*   Open, check, load data from, and close, the initial and final state
*   .csl  files. These files are then merged to one file.
*
      CALL MRGCSL(FULLNAME)
*
*   Open, check, load data from, and close, the merged .csl  file
*
      CALL SETCSLM
*
*   Read mixing coefficients
*
c      CALL READMIX(FULLNAME,INPCI)
*
*   Test mixing coefficients
*
C      IF (NTEST.GT.1000) CALL TESTMIX
*
*   Open, check, load data from, and close the .iso file
*
      CALL SETISO(trim(startdir) // '/isodata')
*
*   Get the remaining information
*
      CALL GETOSD(FULLNAME)
*
*   Open and append a summary of the inputs to the  .sum  file
*
      if (myid .EQ. 0) then
          iii = len_trim(startdir)
          call sys_chdir(trim(startdir),iii,ierr)
          if (ierr.ne.0) then
              print *, "error changing dir!"
              call exit(1)                                            
          endif
          CALL STRSUM(NAME,INPCI,ILBL)
          iii = len_trim(tmpdir)
          call sys_chdir(trim(tmpdir)//'/'//idstring,iii+4,ierr)
          if (ierr.ne.0) then
              print *, "Error changing dir!"
              call exit(1)                                            
          endif
      endif  !myid = 0
*
*   Set up the table of logarithms of factorials
*
      CALL FACTT
*
*   Proceed with the transition calculation
*
      CALL OSCL(NAME,FULLNAME)
*
*   Print completion message
*
      IF (myid .EQ. 0) THEN
          PRINT *
          PRINT *, 'RTRANSITION_MPI: Execution complete.'
      ENDIF
*=======================================================================
*  Execution finished; Statistics output
*=======================================================================
 
      CALL stopmpi2 (myid, nprocs, host, lenhost,
     &                     ncount1, 'RTRANSITION_MPI')                                 
*=======================================================================
*
      STOP
      END
