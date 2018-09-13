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
*   Modified by Gediminas Gaigalas for new spin-angular integration.   *
*                                                                      *
************************************************************************
*
      PROGRAM BIOTR
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
*
*  Total length of scratch space
*
CGG      PARAMETER (LWORK1 = 100000)
      PARAMETER (LWORK1 = 10000000)
*
      PARAMETER (NLMAX = 20)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)

      DIMENSION WORK(LWORK1),NINSHLI(NLMAX),NINSHLF(NLMAX)
      DIMENSION CISHL(NLMAX*NLMAX),S(NNNW,NNNW)
      DIMENSION CFSHL(NLMAX*NLMAX)

      LOGICAL GETYN, YES
      CHARACTER*24 NAME(2)
    
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
*   Debug flags
*
      NTESTG = 000
      NTESTL = 00000
      NTEST = MAX0(NTESTL,NTESTG)
*
      CALL STARTTIME (ncount1, 'BIOTRA3')
CGG      PRINT *
CGG      PRINT *, 'BIOTRA3: Execution starts'
      write(*,*)
      write(*,*) 'RBIOTRANSFORM'
      write(*,*) 'This program transforms the initial and final wave'
      write(*,*) 'functions so that standard tensor albegra can be'
      write(*,*) 'used in evaluation of the transition parameters'  
      write(*,*) 'Input files:  isodata, name1.c, name1.w, name1.(c)m'
      write(*,*) '              name2.c, name2.w, name2.(c)m         '
      print*,'              name1.TB, name2.TB (optional angular files)'
      write(*,*) 'Output files: name1.bw, name1.(c)bm, '
      WRITE(*,*) '              name2.bw, name2.c(bm)'
      write(*,*) '              name1.TB, name2.TB (angular files)'
      write(*,*) 

      PRINT *, 'Default settings?'
      YES = GETYN ()
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
      IF (YES) THEN
         INPCI = 0
      ELSE
         INPCI = 1
      ENDIF
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
      CALL SETISO ('isodata')
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
      CALL FNAME(NAME)
*
*   Open, check, load data from and close the initial state CSL file.
*
      CALL SETCSLB(NAME(1),NCORE1)
*
*   Transfer the data to the initial state COMMON
*           
      CALL TCSL(1)
*
*   Open, check, load data from and close the final state CSL file.
*
      CALL SETCSLB(NAME(2),NCORE2)
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
      CALL LODRWFI (NAME(1),NTESTG)
*
*   Read the the radial orbitals for the initial state
*
      CALL LODRWFF (NAME(2),NTESTG)
*
*   Calculate the radial overlap matrices
*
      WRITE(*,*)
      WRITE(*,*) ' ******************************************'
      WRITE(*,*) '  Overlap matrix before orbital rotations'
      WRITE(*,*) ' *****************************************'
      WRITE(*,*)

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
      WRITE(*,*)
      WRITE(*,*) ' ****************************************'
      WRITE(*,*) '  Overlap matrix after orbital rotations'
      WRITE(*,*) ' ****************************************'
      WRITE(*,*)

      CALL BRKT
*
*  Write the transformed radial functions to file
*
      CALL RADFILE(NAME)
*
*   Obtain one-electron coupling coefficients for the initial state.
*   The coefficients are dumped on files one kappa in turn and
*   thus the different kappa can be manipulated independently.
*   The interface with the transformation part is in the routine mcp
*   
      CALL GENMCP(NAME(1),1,NTESTG,INPCI)
****** added by Yu Zou, Feb.18,2000 ***********
        DO K = 1,KAMAX
          close(UNIT=80+K,STATUS="DELETE")
        ENDDO
****** added by Yu Zou, Feb.18,2000 ***********
*
*   Obtain one-electron coupling coefficients for the final state.
*   The coefficients are dumped on files one kappa in turn and
*   thus the different kappa can be manipulated independently.
*   The interface with the transformation part is in the routine mcp
*   
      CALL GENMCP(NAME(2),2,NTESTG,INPCI)

        DO K = 1,KAMAX
          close(UNIT=80+K,STATUS="DELETE")
        ENDDO
      CALL STOPTIME (ncount1, 'BIOTRANSFORM')
CGG      PRINT *
CGG      PRINT *, ' BIOTRA3: Execution complete.'

      STOP 
      END
