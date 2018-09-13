		SUBROUTINE INIEST(N,NB,NIV,hmx,jcol,irow,BASIS,
     :                    IBLOCK,JBLOCK)
*-----------------------------------------------------------------------
*  	Routine for providing initial estimates from the diagonal 
*	of the matrix. This way was used by Dvdson in atomic structure 
*       calculations. It should be used to obtain estimates when nothing 
*       else is available.
*
*-----------------------------------------------------------------------
        IMPLICIT REAL*8          (A-H,O-Z)
	dimension basis(*)
	dimension IBLOCK(*)
	pointer (iqap,ap(*)),(iqeig,eigval(*)),(iqvec,vec(*))
	pointer (iqwork,work(*)),(iqiwork,iwork(*)),(iqif,IFAIL(*))
	pointer (iqjsh,jsh(*))
	dimension HMX(*),Irow(*),Jcol(0:*)      

*******************************************************************
  
****** alloc space for 100*100 lower triangular
        NS = min(1000, NB)
!once here        NS = min(800, NB)
        !print*, 'NS=',NS,' NIV=',NIV
	call alloc(iqap,NS*(NS+1)/2,8)
	call dinit(NS*(NS+1)/2,0.d0,ap,1)
	
	call alloc(iqeig,NS,8)
	call alloc(iqvec,NS*NIV,8)
	call alloc(iqwork,8*NS,8)
	call alloc(iqiwork,5*NS,4)
	call alloc(iqif,NS,4)
	call alloc(iqjsh,NS,4)
	icount = 0
	
***** separate upper left block of size NS*NS
	do j = 1,N
	 if (icount .ge. NS) goto 101
	 if (IBLOCK(j) .eq. JBLOCK) then
	   icount = icount + 1
	   jsh(icount) = j
	   istart = jcol(j-1) + 1
	   ish = j - icount 
 100       ap(irow(istart)-ish+(icount-1)*(2*NS-icount)/2)
     :        =hmx(istart)
*           print*, 'ic = ',icount, ' i=', irow(istart)-ish
	   istart = istart+1
* check the block structure for zero elements	   
           if ( istart  .gt. jcol(j)) goto 102
	   do ii = irow(istart-1),irow(istart) 
	      if (IBLOCK(ii) .ne. JBLOCK) ish = ish + 1
	   enddo
	   if ( istart  .gt. jcol(j)) goto 102
	   if ( (irow(istart)-ish) .le. NS ) goto 100
*           goto 100
 102       CONTINUE
	 endif
	enddo
 101	CONTINUE
*	print*, 'calling dspevx'
        abstol = slamch('S')
        CALL DSPEVX('Vectors also','In a range','Lower triangular',
     :          NS,AP,-1.,-1.,1,NIV,abstol,
     :          NFOUND,EIGVAL,VEC,NS,work,iwork,IFAIL,INFO)
           IERR = -ABS(INFO)
           if (ierr .NE. 0) print*, 'iniest ierr =', ierr
*           print '(D14.7,X,I2)', (eigval(i),i, i=1,NIV)

*******************************************************************


*
*       ..Build the Basis. 
*
           CALL DINIT(N*NIV,0.D0,BASIS,1)
           istart = 0
           DO 50 J=1,NIV
*scatter the vectors
             DO IC = 1,NS
	      basis(istart+jsh(ic)) = vec((J-1)*NS+IC)
	     ENDDO
*              call dcopy(NS,vec((J-1)*NS+1),1,BASIS(istart+ish(j)),1) 
              istart = istart + N
 50        CONTINUE
        call dcopy(NIV, EIGVAL,1,BASIS(NIV*N+1),1)
        call dalloc(iqap)
        call dalloc(iqeig)
        call dalloc(iqvec)
        call dalloc(iqwork)
        call dalloc(iqiwork)
        call dalloc(iqif)
        call dalloc(iqjsh)

	RETURN
	END
