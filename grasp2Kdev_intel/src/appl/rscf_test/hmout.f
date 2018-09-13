************************************************************************
*                                                                      *
      SUBROUTINE HMOUT (myid, nprocs, ncf,eav)
*                                                                      *
*   Routine for printing the Hamiltonian matrix.                       *
*                                                                      *
*   Written by Farid A Parpia             Last revision: 21 Dec 1992   *
*   Block Version by Xinghong He          Last revision: 30 Jan 1999   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
*
*      integer*8 nelmnt

      POINTER (PNTEMT,EMT(*))
      POINTER (PIENDC,IENDC(0:*))
      POINTER (PNIROW,IROW(*))
*
      COMMON/HMAT/PNTEMT,PIENDC,PNIROW,NELMNT
*
*
!  This is original, correct, full matrix
!      write(99,*)'++++++++++++++++++++++++++++++++++'
       icount = 1
       ibeg = 1
       do ico = myid + 1, ncf, nprocs
          idiag = iendc(ico)
          do list = ibeg, idiag
             iro = irow(list)
             emtt = emt(list)
             if(list == idiag) emtt = emtt-eav
             write (99,*) icount,iro, emtt
!            write (99,*) 'H(',iro,ico,')= ', emt(list)
             icount = icount + 1
          enddo
          ibeg = idiag + 1
       enddo
!      write(99,*)'++++++++++++++++++++++++++++++++++'
!
!     ibeg = 1
!     do ico = 1, ncf, 4
!        ibeg = iendc(ico - 1) + 1
!        idiag = iendc(ico)
!        do list = ibeg, idiag
!           iro = irow(list)
!           write (96,*) 'H(',iro,ico,')= ', emt(list)
!        enddo
!        !ibeg = idiag + 1
!     enddo

!     ibeg = iendc(1) + 1
!     do ico = 2, ncf, 4
!        ibeg = iendc(ico - 1) + 1
!        idiag = iendc(ico)
!        do list = ibeg, idiag
!           iro = irow(list)
!           write (97,*) 'H(',iro,ico,')= ', emt(list)
!        enddo
!        !ibeg = idiag + 1
!     enddo

!     ibeg = iendc(2) + 1
!     do ico = 3, ncf, 4
!        ibeg = iendc(ico - 1) + 1
!        idiag = iendc(ico)
!        do list = ibeg, idiag
!           iro = irow(list)
!           write (98,*) 'H(',iro,ico,')= ', emt(list)
!        enddo
!        !ibeg = idiag + 1
!     enddo

!      write(99,*) '----------------------------'
!     ibeg = iendc(3) + 1
!     do ico = 4, ncf, 4
!        ibeg = iendc(ico - 1) + 1
!        idiag = iendc(ico)
!        do list = ibeg, idiag
!           iro = irow(list)
!           write (99,*) 'H(',iro,ico,')= ', emt(list)
!        enddo
!        !ibeg = idiag + 1
!     enddo
!      write(99,*) '----------------------------'

      END
