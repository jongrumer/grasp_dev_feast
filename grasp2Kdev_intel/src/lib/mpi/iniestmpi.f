************************************************************************
      SUBROUTINE iniestmpi (nmax, ncf, NIV, BASIS, hmx, jcol, irow)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

*  	Routine for providing initial estimates from the diagonal 
*     of the matrix. This way was used by Dvdson in atomic structure 
*     calculations. It should be used to obtain estimates when nothing 
*     else is available.
*
*       nmax is typically 1000
*
*   MPI version by Xinghong He            Last revision: 29 Jul 1998
*
************************************************************************

      DIMENSION basis(*)
      POINTER (iqap,ap(*)),(iqeig,eigval(*)),(iqvec,vec(*))
      POINTER (iqwork,work(*)),(iqiwork,iwork(*)),(iqif,IFAIL(*))
      DIMENSION HMX(*),Irow(*),Jcol(0:*)

      !INCLUDE 'mpif.h'
      COMMON /mpi/ myid, nprocs, ierr
*-----------------------------------------------------------------------
      NS = MIN (nmax, ncf)

      CALL alloc (iqap, (NS*(NS+1))/2, 8)
      CALL dinit (NS*((NS+1))/2, 0.d0, ap, 1)
	
*  Expand the sparse form to normal form for upper-left sub-matrix
*  Offset for the sparse hmx is ACCUMULATED (COUNTED).
*  Offset for the normal ap  is COMPUTED

      joffspar = 0                        ! offset for sparse form
      DO j = myid + 1, ns, nprocs
         joffnorm = (j * (j-1)) / 2       ! offset for normal form
         DO ir = joffspar + 1, jcol(j)
            ap(irow(ir) + joffnorm) = hmx(ir)
         ENDDO
         ! joffspar = joffspar + jcol(j)     ! For rcimpivu ??? No
         joffspar = jcol(j)
      ENDDO

*  Merge ap from all nodes and then send to all nodes

      CALL gdsummpi (ap, (NS*(NS+1))/2)

      CALL alloc (iqeig,   NS,     8)
      CALL alloc (iqvec,   NS*NIV, 8)
      CALL alloc (iqwork,  8*NS,   8)
      CALL alloc (iqiwork, 5*NS,   4)
      CALL alloc (iqif,     NS,    4)

      ABSTOL = 2*SLAMCH('S')
      CALL DSPEVX ('Vectors also','In a range','Upper triangular',
     :          NS,AP,-1.,-1.,1,NIV,ABSTOL,
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

  302 FORMAT (/'Configuration mixing coefficients:')
  303 FORMAT (1X,1P,6D12.4)


      RETURN
      END
