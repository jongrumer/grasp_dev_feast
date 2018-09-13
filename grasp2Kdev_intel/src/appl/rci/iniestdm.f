*************************************************************************
*
	SUBROUTINE INIESTdm (nmax, ncf, NIV, BASIS, hmx)

*  	Routine for providing initial estimates from the diagonal 
*     of the matrix. This way was used by Dvdson in atomic structure 
*     calculations. It should be used to obtain estimates when nothing 
*     else is available.
*
*   Block version by Xinghong He            Last revision: 18 Jun 1998
*
*************************************************************************
      IMPLICIT REAL*8          (A-H, O-Z)
      DIMENSION basis(*)
      POINTER (iqap,ap(*)),(iqeig,eigval(*)),(iqvec,vec(*))
      POINTER (iqwork,work(*)),(iqiwork,iwork(*)),(iqif,IFAIL(*))
      DIMENSION HMX(*)
      POINTER (PNEVAL, EVALDUM)
      COMMON/EIGVAL/EAV,PNEVAL

*-----------------------------------------------------------------------

      myid = 0
      nprocs = 1

      NS = MIN (nmax, ncf)

      CALL alloc (iqap, NS*(NS+1)/2, 8)
      CALL dinit (NS*(NS+1)/2, 0.d0, ap, 1)
	
*  Get the upper left sub-matrix

      DO i = 1, (ns*(ns+1))/2
         ap(i) = hmx(i)
      ENDDO

      CALL alloc (iqeig,   NS,     8)
      CALL alloc (iqvec,   NS*NIV, 8)
      CALL alloc (iqwork,  8*NS,   8)
      CALL alloc (iqiwork, 5*NS,   4)
      CALL alloc (iqif,     NS,    4)

      CALL DSPEVX ('Vectors also','In a range','Upper triangular',
     :          NS,AP,-1.,-1.,1,NIV,0.d0,
     :          NFOUND,EIGVAL,VEC,NS,work,iwork,IFAIL,INFO)
           IERR = -ABS (INFO)

*  Build the Basis. 

      CALL DINIT (ncf*NIV, 0.D0, BASIS, 1)

*  scatter the vectors

      DO J = 1, NIV
	      CALL dcopy (ns, vec(ns*(j-1)+1),1, basis(ncf*(j-1)+1), 1)
      ENDDO

      CALL dcopy (NIV, EIGVAL, 1, BASIS(NIV*ncf+1), 1)

      CALL dalloc (iqap)
      CALL dalloc (iqeig)
      CALL dalloc (iqvec)
      CALL dalloc (iqwork)
      CALL dalloc (iqiwork)
      CALL dalloc (iqif)

      RETURN
      END
