************************************************************************
*                                                                      *
*                                                                      *
      SUBROUTINE CITRAG(CIIN,NCSF,NCIV,L,NSHL,T,NIN,NTESTG,CIOUT, SCR)
*                                                                      *
*   Calculate the action of the  operator                              *
*                                                                      *
*    T(L) = T(nl) T(n-1 L ) .... T(1 l)                                *
*                                                                      *
*   where                                                              *
*                                                                      *
*    T(i l ) = ( sum(k=0,2l+1) t(i l ) ** k ) (t_ilil ** E(il,il))     *
*    t(i l ) = sum(j.lt.i) t(j,i) E(lj,li)                             *
*    E(lj,li) total symmetric shell excitation operator                *
*                                                                      *
*   Modification of the CITRA routine for MCHF                         *
*                                                                      *
*   =====                                                              *
*   Input                                                              *
*   =====                                                              *
*   CIIN      : Input CI vectors ( Destroyed in the process )          *
*   NCSCF     : Length of CI expansion                                 *
*   NCIV      : Number of CI vectors                                   *
*   L         : L value of shells in excitations                       *
*   NSHL      : Number of shells with this L                           *
*   T         : List of coefficients                                   *
*   NIN       : Number of inactive shells                              *
*                                                                      *
*   ======                                                             *
*   Output                                                             *
*   ======                                                             *
*                                                                      *
*   CIOUT : List of output CI vectors                                  *
*                                                                      *
************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*                                                                       
*  Total length of scratch space                                        
*
      !PARAMETER (LWORK2 = 100000)
*
      PARAMETER (NLMAX = 40)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)

      DIMENSION SCR(NCSF,NCIV)
      DIMENSION SCRtmp(NCSF,NCIV)
      DIMENSION CIIN(NCSF,NCIV),T(NSHL,NSHL)
      DIMENSION CIOUT(NCSF,NCIV)
      DIMENSION CIOUTtmp(NCSF,NCIV)
*
      COMMON/SBDAT1/NSHLP(NLMAX,NLMAX),NSHLPP(NLMAX,NNNW)
     :      /ORB4/NP(NNNW),NAK(NNNW)

      include 'mpif.h'
      integer  myid, nprocs, ierr
      COMMON /mpi/ myid, nprocs, ierr                                             

      NTESTL = 00
      NTEST = MAX(NTESTG,NTESTL)
      NTEST = 00

*
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*) ' ***************'
        WRITE(6,*) ' Entering CITRAG'
        WRITE(6,*) ' ***************'
        WRITE(6,*)
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input CI vectors '
        CALL WRTMAT(CIIN,NCSF,NCIV,NCSF,NCIV)
        WRITE(6,*) ' Transformation matrix T'
        CALL WRTMAT(T,NSHL,NSHL,NSHL,NSHL)
      END IF
*
*. Factor from inactive shells
*
      IF(NIN .NE. 0 ) THEN
        FACTOR = 1.0D0
        DO IIN = 1, NIN
          FACTOR = FACTOR*T(IIN,IIN)
        ENDDO
*
*       IPOT = 2*(2*L+1)  (number of m_lm_s. This should be replaced
*                          by (2j+1) corresponding to L)
*
        IPOT = 2*IABS(NAK(NSHLP(L,IIN)))
        FACTOR = FACTOR ** IPOT
        CALL SCALVE(CIIN,FACTOR,NCIV*NCSF)
      END IF
      IF(NIN.EQ.NSHL) CALL COPVEC(CIIN,CIOUT,NCIV*NCSF)
*
      DO I = NIN+1,NSHL
        IF(NTEST.GE.100) WRITE(6,*)  ' Loop I,L = ', I,L
*
*. The diagonal contribution
*
        TII  = T(I,I)
*R Note that the global summation is performed INSIDE TIINIG
        CALL TIINIG(CIIN,NCSF,NCIV,I,L,TII,CIOUT,NTESTG)
c        IF (LWORK2.LT.NCIV*NCSF) THEN
c          WRITE(*,*) 'In CITRAG: Dimension of LWORK2 must be',
c     &               'increased to at least',NCIV*NCSF
c        ENDIF
        CALL COPVEC(CIOUT,SCR,NCIV*NCSF)
*
*.  Off diagonal contributions
*
        XNFACI = 1.0D0
        DO N= 1,2*IABS(NAK(NSHLP(L,I)))
          IF(NTEST.GE.100) WRITE(6,*) ' Loop N = ', N
*
*   T ** (N-1) is supposed to be in SCR, copy to CIIN
*   and apply S
*
          CALL COPVEC(SCR,CIIN,NCIV*NCSF)
          CALL TI1TV(CIIN,NCSF,NCIV,I,L,T(1,I),NSHL,SCRtmp,NTESTG)
*R Remove next statement
*R      CALL COPVEC(SCRtmp,SCR,NCIV*NCSF)
*R Remove pervious statement
          CALL MPI_ALLREDUCE(SCRtmp(1,1),SCR(1,1),NCIV*NCSF,
     &      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)                
          XNFACI = XNFACI/FLOAT(N)
          CALL VECSUM(CIOUT,CIOUT,SCR,1.0D0,XNFACI,NCIV*NCSF)
        ENDDO
        CALL COPVEC(CIOUT,CIIN,NCIV*NCSF)
*
      ENDDO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output CI vectors L = ',L
        CALL WRTMAT(CIOUT,NCSF,NCIV,NCSF,NCIV)
      END IF
*
      RETURN
      END
