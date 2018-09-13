************************************************************************
*                                                                      *
      SUBROUTINE HMOUT (myid, nprocs, ncf)
*                                                                      *
*   Routine for printing the Hamiltonian matrix.                       *
*
*   In this MPI version, each node will print its own part of the 
*   matrix.
*                                                                      *
*   Written by Farid A Parpia             Last revision: 21 Dec 1992   *
*   MPI Version by Xinghong He            Last revision: 30 Jan 1999   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
*
      POINTER (PNTEMT,EMT(*))
      POINTER (PIENDC,IENDC(0:*))
      POINTER (PNIROW,IROW(*))
*
      COMMON/HMAT/PNTEMT,PIENDC,PNIROW,NELMNT
*
*
      ibeg = 1
      do ico = myid + 1, ncf, nprocs
         idiag = iendc(ico)
         do list = ibeg, idiag
            iro = irow(list)
            write (99,*) 'H(',iro,ico,')= ', emt(list)
         enddo
         ibeg = idiag + 1
      enddo

      END
