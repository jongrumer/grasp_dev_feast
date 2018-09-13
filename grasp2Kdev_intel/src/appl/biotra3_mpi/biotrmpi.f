************************************************************************
*                                                                      *
*     This program transforms the initial and final state radial       *
*     orbitals to a biorthonormal set. The CI coefficients are         *
*     then counter transformed as to leave both the initial and        *
*     final state wave functions invariant                             *
*                                                                      *
*     General references:                                              *
*                                                                      *
*     P.A. Malmqvist, Int.J. of Quantum Chemistry, XXX, 479-94 (1986)  *
*     J. Olsen, M.R. Godefroid, P. Jonsson, P.A. Malmqvist and         *
*     C. Froese Fischer, Phys. Rev. E, 4499 (1995)                     *
*                                                                      *
*                                                                      *
*     Program written by                                               *
*                                                                      *
*     Per Jonsson, Department of Computer Science                      *
*                  Vanderbilt University, USA                          *
*                                                                      *
*     Date June 1996                                                   *
*                                                                      *
************************************************************************
*
      PROGRAM BIOTR
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
*
*  Total length of scratch space
*
      PARAMETER (LWORK1 = 100000)
*
      PARAMETER (NLMAX = 40)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590) 
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)

      DIMENSION WORK(LWORK1),NINSHLI(NLMAX),NINSHLF(NLMAX)
      DIMENSION CISHL(NLMAX*NLMAX),S(NNNW,NNNW)
      DIMENSION CFSHL(NLMAX*NLMAX)

      LOGICAL GETYN, YES
      CHARACTER*24 NAME(2)
      CHARACTER*128 FULLNAME(2)
    
      POINTER (PNTRPFII,PFII(NNNP,*)),(PNTRQFII,QFII(NNNP,*))

      POINTER (PNTRPFFF,PFFF(NNNP,*)),(PNTRQFFF,QFFF(NNNP,*))

      COMMON/DEFAULT/NDEF,NDUMP
      COMMON/SBDAT/NAKINVII(NNNW),NSHLII(NLMAX),NSHLPII(NLMAX,NLMAX),   
     :             NAKINVFF(NNNW),NSHLFF(NLMAX),NSHLPFF(NLMAX,NLMAX),   
     :             NSHLPPII(NLMAX,NNNW),NSHLPPFF(NLMAX,NNNW),
     :             NINII(NLMAX),NINFF(NLMAX),IKAPPA(NLMAX),KAMAX
*
*  Initial state commons
*
*     INCLUDE 'initial.inc'
       COMMON/ORB2II/NCFII,NWII
     :      /WAVEII/PZII(NNNW),PNTRPFII,PNTRQFII,MFII(NNNW)
     :      /CIIMAT/CICI(20*NLMAX*NLMAX)
*
*  Final state commons
*
*     INCLUDE 'final.inc'
      COMMON/ORB2FF/NCFFF,NWFF
     :      /WAVEFF/PZFF(NNNW),PNTRPFFF,PNTRQFFF,MFFF(NNNW)
     :      /CFFMAT/CFCI(20*NLMAX*NLMAX)
*
* MPI setup
*

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
     &                  startdir, permdir, tmpdir, 'RBIOTRANSFORM_MPI')
      WRITE (idstring, '(I3.3)') myid
      print*, tmpdir , ' = tmpdir'
 
      lenperm = LEN_TRIM (permdir)
      lentmp = LEN_TRIM (tmpdir)
 
*=======================================================================          
*
*   Debug flags
*
      NTESTG = 000
      NTESTL = 00000
      NTEST = MAX0(NTESTL,NTESTG)
*
      if (myid .eq. 0) then
      PRINT *
*      PRINT *, 'BIOTRA2: Execution starts'
      PRINT * 
      PRINT *, 'Default settings?'
      YES = GETYN ()
      PRINT * 
      IF (YES) THEN
         NDEF = 0
         NDUMP = 1
      ELSE
         NDEF = 1
         PRINT *, 'Dump angular data on file?'
         YES = GETYN ()
         PRINT *
         IF (YES) THEN
            NDUMP = 1
         ELSE
            NDUMP = 0
         ENDIF
      ENDIF
      PRINT *, 'Input from a CI calculation?'
      YES = GETYN ()
      PRINT *
      IF (YES) THEN
         INPCI = 0
      ELSE
         INPCI = 1
      ENDIF
      endif !myid=0
      CALL MPI_Bcast (NDEF,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (NDUMP,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (INPCI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)  
*     
*   Perform machine- and installation-dependent setup
*
      CALL SETMC
*
*   Set up the physical constants
*
      CALL SETCON
*
*   Open, check, load data from, and close the  .iso  file
*
      CALL SETISO (trim(startdir)//'/isodata')
*
*   Determine the parameters controlling the radial grid
*   
      CALL RADPAR
*
*   Generate the radial grid
*
      CALL RADGRD
*
*   Set up the coefficients for the numerical procedures
*
      CALL SETQIC
*
*   Obtain the names of the initial and final state files
*   and open files where the transformed orbitals and CI
*   coefficients are to be dumped
*
      if (myid .eq. 0) CALL FNAME(NAME)
      CALL MPI_Bcast (NAME,48,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      print *, 'DBG: myid=',myid, ' name(1)=', name(1)
      print *, 'DBG: myid=',myid, ' name(2)=', name(2)
      FULLNAME(1) = TRIM(startdir) // '/' // TRIM(NAME(1))
      FULLNAME(2) = TRIM(startdir) // '/' // TRIM(NAME(2))
      print *, 'DBG: myid=',myid, ' startdir =', startdir
      print *, 'DBG: myid=',myid, ' fullname(1)=', fullname(1)
      print *, 'DBG: myid=',myid, ' fullname(2)=', fullname(2)                    
*
*   Open, check, load data from and close the initial state CSL file.
*
      CALL SETCSLB(FULLNAME(1),NCORE1)
*
*   Transfer the data to the initial state COMMON
*           
      CALL TCSL(1)
*
*   Open, check, load data from and close the final state CSL file.
*
      CALL SETCSLB(FULLNAME(2),NCORE2)
*
*   Transfer the data to the final state COMMON
*           
      CALL TCSL(2)
*
*   Determine the number of kappa quantum numbers and
*   the number of orbitals for each kappa quantum number 
*   for the initial state and final states
*
      CALL KAPDATA (NTESTG,NCORE1,NCORE2)
*
*   Read the the radial orbitals for the initial state
*
      CALL LODRWFI (FULLNAME(1),NTESTG)
*
*   Read the the radial orbitals for the initial state
*
      CALL LODRWFF (FULLNAME(2),NTESTG)
*
*   Calculate the radial overlap matrices
*
      if (myid .eq. 0) then
      WRITE(*,*)
      WRITE(*,*) ' ******************************************'
      WRITE(*,*) '  Overlap matrix before orbital rotations'
      WRITE(*,*) ' *****************************************'
      WRITE(*,*)
      endif!myid=0

      CALL BRKT 

      CALL GETS(S,NWII,NWFF)

*
*   Once we have the overlap matrices
*   we can manipulate the initial and final state separately.
*
      MXL = KAMAX
*
*. Calculate biorthonormal orbitals, and orbital matrix
*. for counter transformation of CI coefficients.
*
      CALL BIOTR1(PFII(1,1),QFII(1,1),NSHLII(1),NINSHLI(1),
     :            PFFF(1,1),QFFF(1,1),NSHLFF(1),NINSHLF(1),
     :            NNNP,KAMAX,WORK(1),LWORK1,NTESTG,
     :            CISHL(1),CICI(1),CFSHL(1),CFCI(1) )
      if (myid .eq. 0) then
      WRITE(*,*)
      WRITE(*,*) ' ****************************************'
      WRITE(*,*) '  Overlap matrix after orbital rotations'
      WRITE(*,*) ' ****************************************'
      WRITE(*,*)
      endif !myid=0

      CALL BRKT
*
*  Write the transformed radial functions to file
*
      CALL RADFILE(FULLNAME)
*
*   Obtain one-electron coupling coefficients for the initial state.
*   The coefficients are dumped on files one kappa in turn and
*   thus the different kappa can be manipulated independently.
*   The interface with the transformation part is in the routine mcp
*   
      iii = len_trim(startdir)   
      call sys_chdir(trim(startdir),iii,ierr)
      CALL SETCSLA(NAME(1), ncore_not_used)   
      iii = len_trim(tmpdir)   
      call sys_chdir(trim(tmpdir)//'/'//idstring,iii+4,ierr)
      CALL FACTT
      CALL MCPOUT(NAME(1),startdir,1,NTESTG,INPCI) 
      CALL MCPIN(NAME(1),startdir,1,NTESTG,INPCI) 

****** added by Yu Zou, Feb.18,2000 ***********
*R      DO K = 1,KAMAX
*R        close(UNIT=80+K,STATUS="DELETE")
*R      ENDDO
****** added by Yu Zou, Feb.18,2000 ***********
*
*   Obtain one-electron coupling coefficients for the final state.
*   The coefficients are dumped on files one kappa in turn and
*   thus the different kappa can be manipulated independently.
*   The interface with the transformation part is in the routine mcp
*   
      iii = len_trim(startdir)   
      call sys_chdir(trim(startdir),iii,ierr)
      CALL SETCSLA(NAME(2), ncore_not_used)   
      iii = len_trim(tmpdir)   
      call sys_chdir(trim(tmpdir)//'/'//idstring,iii+4,ierr)
      CALL MCPOUT(NAME(2),startdir,2,NTESTG,INPCI) 
      CALL MCPIN(NAME(2),startdir,2,NTESTG,INPCI) 

*R      DO K = 1,KAMAX
*R        close(UNIT=80+K,STATUS="DELETE")
*R      ENDDO

*=======================================================================
*  Execution finished; Statistics output
*=======================================================================
 
      CALL stopmpi2 (myid, nprocs, host, lenhost,
     &                     ncount1, 'RBIOTRANSFORM_MPI')
*=======================================================================
*                                                                                 

      STOP 
      END

