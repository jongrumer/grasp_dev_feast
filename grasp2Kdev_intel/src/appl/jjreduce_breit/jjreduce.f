************************************************************************
*                                                                      *
      PROGRAM JJREDUCE
*                                                                      *
*   From a set of CSLs this program identifies the ones that           *
*   interact with a given multireference                               *
*                                                                      *
*   This program is a slight modification of the GENMCP program        *
*                                                                      *
*   Written by Per Jonsson                                             *
*                                                                      *
*   Modified by Gediminas Gaigalas for new spin-angular integration.   *
*                                                                      *
************************************************************************
*
      LOGICAL RESTRT,DEBUG,GETYN,YES
      COMMON/DEFAULT/NDEF
      INTEGER HAMILTONIAN
      CHARACTER(LEN=1500) STRING
*
      write(*,*)
      write(*,*) 'RCSFINTERACT'
      write(*,*) 'This program determines all the CSFs in rcsf.inp that'
      write(*,*) 'interact with CSFs in the multireference rcsfmr.inp'
      write(*,*) 'Interacting CSFs written to rcsf.out'
      write(*,*) 'Input files: rcsf.inp, rcsfmr.inp'
      write(*,*) 'Output file: rcsf.out'
      write(*,*)

      PRINT *, 'Default settings?'
      YES = GETYN ()
      IF (YES) THEN
         NDEF = 0
      ELSE
         NDEF = 1
      ENDIF

      WRITE(*,*) 'Reduction based on Dirac-Coulomb (1) or'
      WRITE(*,*) 'Dirac-Coulomb-Breit (2) Hamiltonian?'
      READ(*,*) HAMILTONIAN
      WRITE(*,*) 
*
*     The rcsf.inp and rcsfmr.inp file are in block form.
*     Use readwrite and a scratch file to create rcsf.inp and rcsfmr.inp
*     in non-block form
*

      open(unit=195,file='rcsf.inp',status='old',form='formatted',
     :          action='readwrite')

      open(unit=196,status='scratch',form='formatted')

      do 
          read(195,'(a)',end=991) string
          write(196,'(a)') trim(string)
      end do

991   continue

      rewind(195)
      rewind(196)

      do
          read(196,'(a)',end=992) string
          if (string(1:2).ne.' *') write(195,'(a)') trim(string)
      end do

992   continue

      close(195)

      open(unit=197,file='rcsfmr.inp',status='old',form='formatted',
     :          action='readwrite')

      open(unit=198,status='scratch',form='formatted')

      do
          read(197,'(a)',end=993) string
          write(198,'(a)') trim(string)
      end do

993   continue

      rewind(197)
      rewind(198)

      do
          read(198,'(a)',end=994) string
          if (string(1:2).ne.' *') write(197,'(a)') trim(string)
      end do

994   continue

      close(197)

*
*   Check compatibility of plant substitutions
*
      CALL CHKPLT
*
*   Determine if there is to be any debug printout; this will be
*   made on the  .dbg  file
*
      CALL SETDBG (DEBUG)
*
*   Perform machine- and installation-dependent setup
*
      CALL SETMC
*
*   Open the  .sum  file
*
C      IF (NDEF.NE.0) CALL SETSUM
*
*   Cetermine the names of the file to be reduce and the mr file
*
*   Open, check, load data from, and close the  csl  file
*
      CALL SETCSLB
*
*   Determine where the multireference CSFs reside in the csl file
*
      CALL IDENTY
*
*   Set up the  mcp  files; determine if this is a restart
*
C      CALL SETMCP(RESTRT)
      RESTRT = .FALSE.
*
*   Append a summary of the inputs to the  .sum  file
*
C      IF (NDEF.NE.0) CALL STRSUM
*
*   Set up the table of logarithms of factorials for use by
*   angular modules
*
      write(*,*) 'Before FACTT'
      CALL FACTT
      write(*,*) 'After FACTT'
*
*   Proceed with the generation of MCP coefficients
*
      write(*,*) 'Before MCP'
      CALL MCP(RESTRT,HAMILTONIAN)
      write(*,*) 'After MCP'
*
*   All done; close files that are still open
*
      CLOSE (24)
C      IF (DEBUG) CLOSE (99)
*
*
*   Convert rcsf.inp, rcsfmr.inp to block form again

      open(unit=195,file='rcsf.inp',status='old',form='formatted',
     :          action='readwrite')

      rewind(196)

      do
          read(196,'(a)',end=995) string
          write(195,'(a)') trim(string)
      end do

995   continue

      close(195)
      close(196)

      open(unit=197,file='rcsfmr.inp',status='old',form='formatted',
     :          action='readwrite')

      rewind(198)

      do
          read(198,'(a)',end=996) string
          write(197,'(a)') trim(string)
      end do

996   continue

      close(197)
      close(198)

*
*  Finally, put rcsf.out in block form
*
      call rcsfblock


      write(*,*) 'RCSFINTERACT: Execution complete'

      STOP
      END
