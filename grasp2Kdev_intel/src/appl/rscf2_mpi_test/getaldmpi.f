************************************************************************
*
      SUBROUTINE GETALDmpi
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
*
*   Interactively determines the data governing AL problem.
*
*   Quantities to be determined:
*     ncmin, wt(1:ncf), ucf(1:nw), nscf, nsic, orthst
*
*   Call(s) to: [LIB92]: ALLOC, IQ.
*
*   Written by Farid A. Parpia            Last revision: 19 Dec 1992
*   MPI version by Xinghong He            Last revision: 13 Jul 1998
*
************************************************************************

      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)

      POINTER (pccmin,ICCMIN(*))    ! Not used
      POINTER (pweigh,WEIGHT(*))    ! Not used
      
      INTEGER*4 IQADUM
      POINTER (PNTRIQ,IQADUM)
      LOGICAL getyn, lfix, orthst, yes

      POINTER (pntrwt,wt(*))

      COMMON/def4/accy,nscf,nsic,nsolv
     :      /def5/pntrwt,pweigh
     :      /def7/pccmin,ncmin,ncmax
     :      /fixd/nfix,lfix(NNNW)
     :      /orb2/ncf,nw,pntriq
     :      /orthct/orthst
     :      /scf1/ucf(NNNW)

      COMMON/iounit/istdi,istdo,istde

      INCLUDE 'mpif.h'
      COMMON /mpi/ myid, nprocs, ierr
!-----------------------------------------------------------------------
      IF (myid .EQ. 0) THEN
         WRITE (istde,*) '(E)AL type calculation; H(DC) will not be '
     &,              'diagonalised;'
         WRITE (istde,*) 'getald ...'
         WRITE (istde,*) 'ncf=', ncf
      ENDIF

      CALL alloc (pntrwt, ncf, 8)

      IF (myid .EQ. 0) CALL getaldwt (ncf, wt)

      CALL MPI_Bcast (wt(1), ncf, MPI_DOUBLE_PRECISION, 0, 
     &                         MPI_COMM_WORLD, ierr)

      DO j = 1, nw
         sum = 0.D0
         DO i = 1, ncf
            sum = sum + wt(i) * DBLE (IQ (j,i))
         ENDDO
         ucf(j) = sum
      ENDDO

      ncmin = 0
      nscf = 12
      nsic = 2 + (nw - nfix) / 4
      orthst = .FALSE.

      RETURN
      END
