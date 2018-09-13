************************************************************************
*                                                                      *
      SUBROUTINE RADPAR
*                                                                      *
*   This subroutine sets the parameters controlling the radial grid    *
*                                                                      *
*   Last revision: June 1996                                           *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
*
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590) 
CGG      PARAMETER (NNN1 = NNNP+10)

      LOGICAL GETYN,YES
      CHARACTER*1 ANSWER
      COMMON/NPAR/PARM(2),NPARM
     :      /DEF1/EMN,IONCTY,NELEC,Z
     :      /DEF2/C
     :      /DEF4/ACCY,NSCF,NSIC,NSOLV
     :      /DEF9/CVAC,PI
     :      /DEFAULT/NDEF
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N

      INCLUDE 'mpif.h'
      INTEGER myid, nprocs, ierr
      COMMON /mpi/ myid, nprocs, ierr                                             

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
        if (myid .eq. 0) then
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
        endif !myid=0
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

      IF (NDEF.NE.0) THEN
       if(myid .EQ. 0) then
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
       endif !myid=0
         CALL MPI_Bcast (C, 1, MPI_DOUBLE_PRECISION, 0,
     &   MPI_COMM_WORLD, ierr) 
      ELSE
         C = CVAC
      ENDIF
      
      RETURN
      END
