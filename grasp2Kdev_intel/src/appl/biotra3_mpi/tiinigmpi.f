************************************************************************
*                                                                      *
      SUBROUTINE TIINIG(CIIN,NCSF,NCIV,I,L,CONST,CIOUT,NTESTG)
*                                                                      *
*   Calculates the action of the operator                              *
*   Const ** E(li,li) on a set of vectors                              *
*                                                                      *
*   Adapted for GRASP, daughter of TIINI, born February 1996           *
*                                                                      *
*   =====                                                              *
*   Input                                                              *
*   =====                                                              *
*   CIIN      : Input CI vectors                                       *
*   NCSCF     : Length of CI expansion                                 *
*   NCIV      : Number of CI vectors                                   *
*   I         : Shell number                                           *
*   L         : symmetry                                               *
*   NSHLP(L,K): Gives the shell number as defined in getcsl for        *
*               the K:th shell with symmetry L                         *
*   CONST     : The constant                                           *
*                                                                      *
*   ======                                                             *
*   Output                                                             *
*   ======                                                             *
*                                                                      *
*   CIOUT : List of output CI vectors                                  *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'parameters.def'
CGG      PARAMETER (KEY = 121)
CGG begin
      PARAMETER (KEY = KEYORB)
CGG end
      PARAMETER (NLMAX = 40)                                            
CGG      PARAMETER (NNNW = 120)

      DIMENSION CIIN(NCSF,NCIV)
      DIMENSION CIOUT(NCSF,NCIV)
      DIMENSION CIOUTtmp(NCSF,NCIV)

      POINTER (PJANN,JANN(*)),(PPPINT,INTGRL(*)),
     :        (PJBNN,JBNN(1)),(PCNN,CNN(1)),(PINTPTR,INTPTR(1))

      COMMON/MCPDATA/PJANN,PJBNN,PPPINT,PCNN,PINTPTR,NCOEFF,NINT
      COMMON/SBDAT1/NSHLP(NLMAX,NLMAX),NSHLPP(NLMAX,NNNW)
      include 'mpif.h'
      integer  myid, nprocs, ierr
      COMMON /mpi/ myid, nprocs, ierr                                             
*
      NTESTL = 00
      NTEST = MAX(NTESTL,NTESTG)


      IF(NTEST.GE.10) WRITE(6,*) ' Entering TIINI'

      CALL SETVEC(CIOUT,0.0D0,NCSF*NCIV)
      if (NCOEFF .EQ. 0) goto 200
*
*.  Obtain address of first coupling coefficient for h(il,il) :  IFIRST
*.  Obtain number of        coupling coefficient for h(il,il) :  NFOUND
*.
*     NFOUND : Number of coefficients obtained
*     IVAL   : actual RACAH coefficient <CSF(L)!E(il,il)!CSF(R)>
*     ILEFT  = CSF(L) ?
*
      NFOUND = 0
      DO K = 1,NINT
        IA = INTGRL(K)/KEY
        IB = MOD(INTGRL(K),KEY)
        IF (NSHLP(L,I).EQ.IA.AND.IA.EQ.IB) THEN
          IF (K.EQ.1) THEN
            IFIRST = 1
          ELSE
            IFIRST = INTPTR(K-1) + 1
          ENDIF
          NFOUND = INTPTR(K) - IFIRST + 1
          GOTO 15
        ENDIF
      ENDDO

   15 DO IELMNT = 1, NFOUND
        IVAL = CNN(IFIRST-1+IELMNT)
        CONSTN = CONST ** IVAL
        ILEFT = JANN(IFIRST-1+IELMNT)
        DO IVEC=1, NCIV
          CIOUT(ILEFT,IVEC) = CONSTN * CIIN(ILEFT,IVEC)
        ENDDO
      ENDDO

 200  call copvec(CIOUT, CIOUTtmp, NCIV*NCSF)
      call MPI_ALLREDUCE(CIOUTtmp(1,1), CIOUT(1,1), NCIV*NCSF,
     &  MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
*
*.  The previous provided us with all
*   terms with nonvanishing occupation.
*   For terms with vanishing occupation of il,
*   just copy coefficients, since (x) ** 0 = 1
*
      DO IVEC = 1, NCIV
        DO IELMNT = 1, NCSF
          IF(CIOUT(IELMNT,IVEC).EQ.0.0D0)
     &     CIOUT(IELMNT,IVEC) = CIIN(IELMNT,IVEC)
        ENDDO
      ENDDO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Input and output vectors from TIINI I,L',I,L
        CALL WRTMAT(CIIN,NCSF,NCIV,NCSF,NCIV)
        WRITE(6,*)
        CALL WRTMAT(CIOUT,NCSF,NCIV,NCSF,NCIV)
        WRITE(6,*)
      END IF
*
      IF(NTEST.GE.10) WRITE(6,*) ' Leaving  TIINI'

      RETURN
      END
