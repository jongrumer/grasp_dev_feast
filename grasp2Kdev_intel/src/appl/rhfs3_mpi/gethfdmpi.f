************************************************************************
*                                                                      *
      SUBROUTINE GETHFDmpi (NAME)
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
c
cbieron include 'parameters.def'
c
c
      include 'parameters.def'
c
c      PARAMETER (NNNP = 390)
c      PARAMETER (NNN1 = 400)
c      PARAMETER (NNNW = 120)
c
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL GETYN,LFORDR,LTRANS,LVP,LSE,LNMS,LSMS,YES
      CHARACTER*24 NAME
cbieron
cb
      CHARACTER*24 name_rwf
*
      POINTER (PNTRPF,PF(NNNP,*)),(PNTRQF,QF(NNNP,*))
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
      INCLUDE 'mpif.h'
      INTEGER myid, nprocs, ierr, lenhost
      COMMON /mpi/ myid, nprocs, ierr


*   Open, check, load data from, and close the  .iso  file
*
cb CALL SETISOmpi moved to hfs92mpi.f
cb     CALL SETISOmpi ('isodata')
*
      if (myid .EQ. 0) then

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
      IF (NDEF.NE.0) THEN
*
         PRINT *, 'Treat contributions of some CSFs'
         PRINT *, ' as first-order perturbations?'
         YES = GETYN ()
         IF (YES) THEN
            LFORDR = .TRUE.
            PRINT *, 'The contribution of CSFs'
            PRINT *, ' 1 -- ICCUT will be treated'
            PRINT *, ' variationally; the remainder'
            PRINT *, ' perturbatively; enter ICCUT:'
            READ *, ICCUT
         ELSE
            LFORDR = .FALSE.
            ICCUT = 0
         ENDIF
      ELSE
         LFORDR = .FALSE.
         ICCUT = 0
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

      endif   !myid=0

      CALL MPI_Bcast (C, 1, MPI_DOUBLE_PRECISION, 0,
     &                          MPI_COMM_WORLD, ierr)
      CALL MPI_Bcast (LFORDR, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_Bcast (ICCUT,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (RNT, 1, MPI_DOUBLE_PRECISION, 0,
     &                          MPI_COMM_WORLD, ierr)
      CALL MPI_Bcast (H, 1, MPI_DOUBLE_PRECISION, 0,
     &                          MPI_COMM_WORLD, ierr)
      CALL MPI_Bcast (N,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (HP, 1, MPI_DOUBLE_PRECISION, 0,
     &                          MPI_COMM_WORLD, ierr)

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
*   Generate $- r \times V_nuc (r)$
*
      CALL NUCPOT
*
*   Load the radial wavefunctions
*
cbieron 
cb SETRWFA moved to hfs92mpi.f
cb     NAME_len = len_trim(NAME)  
cb     name_rwf = NAME(1:NAME_len)//'.w'
cb      WRITE (istde,*) name_rwf
cb     CALL SETRWFA(name_rwf)
*
      RETURN
      END
