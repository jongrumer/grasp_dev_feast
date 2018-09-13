! Note:
! HIEND, ISELEC() not used outside dvdson

      subroutine GDVD( OP, N, LIM, DIAG, ILOW, IHIGH,
     :                 ISELEC, NIV, MBLOCK, CRITE, CRITC, CRITR,
     :                 ORTHO, MAXITER,
     :                 WORK, IWRSZ, IWORK, IIWSZ, 
     :                 HIEND, NLOOPS,NMV,IERR
     :               )

      IMPLICIT REAL*8           (A-H,O-Z)
      logical HIEND

      dimension IREV(7),WORK(IWRSZ),IWORK(IIWSZ),ISELEC(LIM)
      dimension DIAG(N)
      EXTERNAL OP
***********************************************************************
*  NOC = number of orthogonalization constraints
*  IRC = reverse communication switch

      NOC = 0
      IRC = 0

************************************************************************
* CALLING DAVIDSON with reverse communication
************************************************************************
* Initial estimates
* 
      if (NIV .EQ. 0) then
         print *, 'GDVD Error : No initial estimate!!!'
         IERR = -1000
         return 
      endif

      !ttt=etime_(tarray)
      nmv=0
      !xhh print*, 'MBLOCK = ', mblock
      !xhh print *, ' gdvd:  niv = ', niv
 99   CALL DVDSON(IRC,IREV,N,LIM,NOC,
     :             ILOW,IHIGH,ISELEC,NIV,MBLOCK,
     :             CRITE,CRITC,CRITR,MAXITER,
     :              WORK,IWRSZ,IWORK,IIWSZ,
     :             HIEND,NLOOPS,IERR)
*
* * * Start Reverse Communication * * * * * * * * * * * * * * * * * * * * 
*
      NB   = IREV(1)
      iw1  = IREV(2)
      iw2  = IREV(3)
      iw3  = IREV(4)
      iiw  = IREV(5)
      iin  = IREV(6)
      iw4  = IREV(7)
 
      if (IRC.eq.1) then
********** ..Preconditioning. Solve NB times(M work(iw2)=work(iw1))
*          ..Results always on work(iw2) 

         icur=0
         do 100 i=1,NB
            indx   = IWORK(iiw+i-1) - 1
            value  = WORK( iw3 + indx )
            ovalue = WORK( iw4 + indx )
            rnorm  = WORK( iin + indx )
            epsil  = WORK( iin + LIM + indx )
* The current approximation of the eigenvector is x=Bc
* If needed it should be saved e.g: call dcopy(N,work(iw2),1,curx,1)

*            write(*,11) nloops,value,rnorm
  11         format('It ',I4,'  Dl',D10.3, ' Res:',D10.3)
*
*-------------
* (M-lI)      .. Compute temporarily (M-valueI) (M preconditioner)
*             .. Needed for (M-valueI)^-1 res
*             .. Here M is the DIAG
            do j=1,N
               DIAG(j)=DIAG(j)-value
            enddo
*-------------
** Choice of Diagonal preconditioning
            do j=1,N
               if (abs(DIAG(j)).gt.(1D-6)) then
                  WORK(iw2+icur+j-1)=WORK(iw1+icur+j-1)/DIAG(j)
               else 
                  WORK(iw2+icur+j-1)=WORK(iw1+icur+j-1)*1.D6
               endif
            enddo
**-------------
** e.g: For No preconditioner: Lanczos
**            call dcopy(N,WORK(iw1+icur),1,WORK(iw2+icur),1)
**-------------
** (M+lI)      .. Restore (M+valueI)
            do j=1,N
               DIAG(j)=DIAG(j)+value
            enddo

            icur=icur+N
 100     continue

         goto 99

**********
      else if ((IRC.eq.2).OR.(IRC.eq.3)) then
********** ..Matrix-vector multiply. 
         call OP( N,NB,work(iw1),work(iw2) )
         nmv=nmv+NB

         goto  99
      endif
* * * * * End of Reverse Communication * * * * * * * * * * * * * * * * *

      return
      end
