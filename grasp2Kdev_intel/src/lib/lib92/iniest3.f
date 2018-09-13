************************************************************************
*
      SUBROUTINE iniest3 (nmax, ncf, NIV, evmin, evmax, hmx, jcol, irow)

*  Serial version of iniestmpi.
*  Structure of the input sparse matrix hmx:
*    . It's a 1-d array
*    . Length: 1 to jcol(ncf)
*    . Number of non-zero elements for column j is:
*          jcol(j) - jcol(j-1) + 1
*    . Row index for element hmx(i) is irow(i)
*  Xinghong He  98-10-28

************************************************************************
      IMPLICIT REAL*8          (A-H,O-Z)
!     DIMENSION basis(*)
      POINTER (iqap,ap(*)),(iqeig,eigval(*)),(iqvec,vec(*))
      POINTER (iqwork,work(*)),(iqiwork,iwork(*)),(iqif,IFAIL(*))
      DIMENSION HMX(*),Irow(*),Jcol(0:*)

      double precision :: evmin, evmax
*-----------------------------------------------------------------------
      NS = MIN (nmax, ncf)
!      PRINT *, 'iniest2: ns=',ns,' niv=',niV

      CALL alloc (iqap, (NS*(NS+1))/2, 8)
      CALL dinit ((NS*(NS+1))/2, 0.d0, ap, 1)
	
*  Expand the sparse form to normal form for upper-right sub-matrix

      joffspar = 0                        ! offset for sparse form
      DO j = 1, ns 
         joffnorm = (j * (j-1)) / 2       ! offset for normal form
         DO ir = joffspar + 1, jcol(j)
            ap(irow(ir) + joffnorm) = hmx(ir)
         ENDDO
         joffspar = jcol(j)
      ENDDO

*  Merge ap from all nodes and then send to all nodes

      CALL alloc (iqeig,   NS,     8)
      CALL alloc (iqvec,   NS*NIV, 8)
      CALL alloc (iqwork,  8*NS,   8)
      CALL alloc (iqiwork, 5*NS,   4)
      CALL alloc (iqif,     NS,    4)

      CALL DSPEVX ('Vectors also','In a range','Upper triangular',
     :          NS,AP,-1.,-1.,1,NIV,0.d0,
     :          NFOUND,EIGVAL,VEC,NS,work,iwork,IFAIL,INFO)
      IERR = -ABS (INFO)

      evmin = minval(eigval(1:niv))
      evmax = maxval(eigval(1:niv))

*  Build the Basis. 

!     CALL DINIT (ncf*NIV, 0.D0, BASIS, 1)

*  scatter the vectors

!     DO J = 1, NIV
!      CALL dcopy (ns, vec(ns*(j-1)+1),1, basis(ncf*(j-1)+1), 1)
!     ENDDO

!     CALL dcopy (NIV, EIGVAL, 1, BASIS(NIV*ncf+1), 1)

      CALL dalloc (iqap)
      CALL dalloc (iqeig)
      CALL dalloc (iqvec)
      CALL dalloc (iqwork)
      CALL dalloc (iqiwork)
      CALL dalloc (iqif)

      RETURN
      END
