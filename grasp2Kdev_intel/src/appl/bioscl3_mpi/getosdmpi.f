************************************************************************
*                                                                      *
      SUBROUTINE GETOSD(NAME)
*                                                                      *
*   Interactively determines the data governing the transition prob-   *
*   lem.                                                               *
*                                                                      *
*   Written by Per Jonsson                Last revision: June 1996     *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)

      include 'parameters.def'
CFF      PARAMETER (NNNP = 590) 
CFF      PARAMETER (NNN1 = NNNP+10)
CFF      PARAMETER (NNNW = 120)

      LOGICAL LFORDR,LTC,LTRANS,LVP,LSE,LNMS,LSMS,GETYN,YES
      CHARACTER*128 NAME(2)
      CHARACTER*4 CUNITS
      CHARACTER*1 ANSWER
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
     :      /OSC7/LTC(10)
      INCLUDE 'mpif.h'
      INTEGER myid, nprocs, ierr
      COMMON /mpi/ myid, nprocs, ierr                                           
*
*   Determine the physical effects specifications
*
      IF (NDEF.NE.0) THEN
         IF (myid .EQ. 0) THEN
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
         ENDIF !myid=0
         CALL MPI_Bcast (C, 1, MPI_DOUBLE_PRECISION, 0,
     &   MPI_COMM_WORLD, ierr)                                
      ELSE
         C = CVAC
      ENDIF

*
      LFORDR = .FALSE.
      ICCUT = 0
*
*   Determine the multipolarity and parity of the transitions
*
      CALL GETRMP
*
*   Determine the units for the printout and other options
*
      DO I = 1,10
        LTC(I) = .FALSE.
      ENDDO
*
      IF (NDEF.NE.0) THEN
       IF (myid .EQ.0) THEN
        PRINT *, 'Which units are to be used to'
        PRINT *, ' express the transition energies?'
        PRINT *, '    A    : Angstrom:'
        PRINT *, '    eV   : electron volts;'
        PRINT *, '    Hart : Hartree atomic units;'
        PRINT *, '    Hz   : Hertz;'
        PRINT *, '    Kays : Kaysers [cm**(-1)];'
    2   READ (*,'(A)') CUNITS
        IF (CUNITS(1:1) .EQ. 'A') THEN
          LTC(1) = .TRUE.
        ELSEIF (CUNITS(1:2) .EQ. 'eV') THEN
          LTC(2) = .TRUE.
        ELSEIF (CUNITS(1:4) .EQ. 'Hart') THEN
          LTC(3) = .TRUE.
        ELSEIF (CUNITS(1:2) .EQ. 'Hz') THEN
          LTC(4) = .TRUE.
        ELSEIF (CUNITS(1:4) .EQ. 'Kays') THEN
          LTC(5) = .TRUE.
        ELSE
          PRINT *, 'GETOSD: Unable to interpret string;'
          PRINT *, ' reenter ...'
          GOTO 2
        ENDIF
       ENDIF !myid=0
       CALL MPI_Bcast(LTC(1),5,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      ELSE
        LTC(5) = .TRUE.
      ENDIF
*
      !IF (myid .EQ. 0) THEN
      !    PRINT *, 'Sort transitions by energy?'
      !    YES = GETYN ()
      !    IF (YES) LTC(6) = .TRUE.
      !ENDIF !myid=0
      CALL MPI_Bcast(LTC(6),1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      
*
      IF (NDEF.NE.0) THEN
       IF (myid .EQ. 0) THEN
        PRINT *, 'Einstein A and B coefficients are'
        PRINT *, ' printed in SI units; use Hartree'
        PRINT *, ' atomic units instead?'
        YES = GETYN ()
        IF (YES) LTC(7) = .TRUE.
       ENDIF !myid=0
       CALL MPI_Bcast(LTC(7),1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      ELSE
        LTC(7) = .FALSE.
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
        IF (myid .EQ. 0) THEN
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
        ENDIF !myid=0
        CALL MPI_Bcast (RNT, 1, MPI_DOUBLE_PRECISION,
     &      0,MPI_COMM_WORLD,ierr) 
        CALL MPI_Bcast (H, 1, MPI_DOUBLE_PRECISION,
     &      0,MPI_COMM_WORLD,ierr) 
        CALL MPI_Bcast (HP, 1, MPI_DOUBLE_PRECISION,
     &      0,MPI_COMM_WORLD,ierr) 
        CALL MPI_Bcast (N, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr) 
      ENDIF
*
*   ACCY is an estimate of the accuracy of the numerical procedures
*
      ACCY = H**6
*
*   Set up the coefficients for the numerical proceduRes
*
      CALL SETQIC
*
*   Generate the radial grid and all associated arrays
*
      CALL RADGRD
*
*   Load the initial state radial wavefunctions. 
*
      CALL LODRWFI(NAME(1))
*
*   Load the final state radial wavefunctions. 
*
      CALL LODRWFF(NAME(2))

*   Construct the radial overlap matrix. 
*
      CALL BRKT

      RETURN
      END
