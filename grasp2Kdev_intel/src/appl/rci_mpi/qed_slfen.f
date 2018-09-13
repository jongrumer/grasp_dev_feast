************************************************************************
*                                                                      *
      SUBROUTINE QED_SLFEN (SLFINT)
*                                                                      *
*   This  routine estimates the F(Z\alpha) function of self energy for *
*   each orbital.                                                      *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, DALLOC, IQ, QUAD, RALLOC, SHIELD.      *
*               [RCI92]: FZALF, HOVLAP.                                *
*                                                                      *
*   Modified from subroutine QED by Yu Zou, Last update: 13 Mar 2000   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
      POINTER (PNTRIQ,RIQDUMMY)
      POINTER (PNIVEC,IVECDUMMY)
      CHARACTER*2 NH, npchar*1, nakchar
*
*     DIMENSION UCF(*)
*
      POINTER (PNEVEC,EVEC(*))
      POINTER (PNTRPF,PF(NNNP,*))
      POINTER (PNTRQF,QF(NNNP,*))
*
      INCLUDE 'mpif.h'
      INTEGER myid, nprocs, ierr, lenhost
      COMMON /mpi/ myid, nprocs, ierr
*
      DIMENSION PTEMP(NNNP),QTEMP(NNNP),SLFINT(NNNW)
*
      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /DEF2/C
     :      /DEF9/CVAC,PI
     :      /EIGVEC/PNEVEC
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /HORB/PH(NNNP),QH(NNNP)
     :      /NPAR/PARM(2),NPARM
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /PRNT/NVEC,PNIVEC,NVECMX
     :      /TATB/TA(NNN1),TB(NNN1),MTP
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
     :      /DEFAULT/NDEF
     :      /QEDCUT/NQEDCUT,NQEDMAX
!
! Pre-set tolerable number for iteration in finding effective
! nuclear charge. 
!
      maxiter = 20
!
      IF (myid .EQ. 0) THEN
C       IF (NDEF .NE. 0) THEN
C         WRITE(*,*) 'In QED_SLFEN: set maximal principal quantum number'
C         WRITE(*,*) 'for including self-energy for orbital. Number must'
C         WRITE(*,*) 'be 8 or less'
C         READ(*,*) NPJMAX
C       ELSE
C         NPJMAX = 8
C       END IF
        IF (NQEDCUT.EQ.1) THEN
           NPJMAX = NQEDMAX
        ELSE
           NPJMAX = 8
        END IF
      END IF
      CALL MPI_Bcast (NPJMAX,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!
*
      DO 14 J = 1,NW
*
         NPJ = NP(J)
*
         IF (NPJ .LE. NPJMAX) THEN
*
*   Only orbitals with principal quantum number 8 or less can
*   be treated by this section of code
*
            KAPPA = NAK(J)
*
*   Begin by transferring the function to a temporary array
*
            MFJ = MF(J)
*
            PTEMP(1) = 0.0D 00
            QTEMP(1) = 0.0D 00
            DO 5 I = 2,MFJ
               PTEMP(I) = PF(I,J)
               QTEMP(I) = QF(I,J)
    5       CONTINUE
            zeff =z
            ratio=ratden(PTEMP,QTEMP,MFJ,NPJ,KAPPA,ZEFF)
            VALU =  ratio
     :             *FZALF (NPJ,KAPPA,ZEFF)
     :             /DBLE (NPJ**3)
   13       SLFINT(J) = VALU*ZEFF**4/(PI*C**3)
*
         ELSE
*
*   The self-energy for orbitals with principal quantum number
*   greater than 8 is set to zero
*
            SLFINT(J) = 0.0D 00
*
         ENDIF
*
   14 CONTINUE
*
*   Deallocate storage for the `generalised occupation numbers'
*
*
      RETURN
      END
