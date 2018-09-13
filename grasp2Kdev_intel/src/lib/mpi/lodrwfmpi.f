************************************************************************
      SUBROUTINE LODRWFmpi (ierror)
      IMPLICIT REAL*8          (A-H, O-Z)

*   This subroutine loads radial wavefunctions from the .rwf  file
*   and performs some related setup. It does not handle any error
*   so the caller has to do it.
*       ierror=0, normal
*       ierror=1, error
*
*   Used by rcimpivu, rscfmpi, rcimpi
*
*   Call(s) to: [LIB92]: ALLOC, DALLOC, INTRPQ, ORTHSC.                *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 05 Oct 1992   *
*   MPI version by Xinghong He            Last revision: 27 May 1997   *
*
************************************************************************

CGG      PARAMETER (NNNP = 590)
      PARAMETER (NNNP = 5000)
      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
      PARAMETER (NNNW = 214)

      POINTER (PNTRIQ,RIQDUM)
      LOGICAL LDBPR
      CHARACTER*2 NH
*
      POINTER (PNTRPA,PA(*))
      POINTER (PNTRQA,QA(*))
      POINTER (PNTRRA,RA(*))
*
      POINTER (PNTRPF,PF(NNNP,*))
      POINTER (PNTRQF,QF(NNNP,*))
*
      COMMON/DEBUGR/LDBPR(30)
     :      /DEF1/EMN,IONCTY,NELEC,Z
     :      /DEF2/C
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /NPAR/PARM(2),NPARM
     :      /ORB1/E(NNNW),GAMA(NNNW)
     :      /ORB2/NCF,NW,PNTRIQ
      COMMON/ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)

      INCLUDE 'mpif.h'
      COMMON /mpi/ myid, nprocs, ierr
!-----------------------------------------------------------------------
*
*   Allocate storage to orbital arrays
*
      CALL ALLOC (PNTRPF,NNNP*NW,8)
      CALL ALLOC (PNTRQF,NNNP*NW,8)
*
*   Setup: (1) Orbital arrays to zero
*          (2) Array E to -1 (no orbitals estimated)
*          (3) Parameters GAMMA for each orbital
*
      CON = Z / C
      CON = CON * CON

      DO J = 1, NW

         DO I = 1, N
            PF(I,J) = 0.D0
            QF(I,J) = 0.D0
         ENDDO

         E(J) = -1.D0
         K = ABS (NAK(J))
         IF (NPARM .GT. 0) THEN
            GAMA(J) = DBLE (K)
         ELSEIF (NPARM .EQ. 0) THEN
            FKK = DBLE (K*K)
            IF (FKK .GE. CON) THEN
               GAMA(J) = SQRT (FKK-CON)
            ELSE
               !WRITE (istde,*) 'LODRWF: Imaginary gamma parameter'
               !WRITE (istde,*) ' for ',NP(J),NH(J),' orbital; the'
               !WRITE (istde,*) ' point model for the nucleus'
               !WRITE (istde,*) ' is inappropriate for Z > ',C,'.'
               CALL stopmpi ('lodrwfmpi: Inappropriate gamma', myid)
            ENDIF
         ENDIF
      ENDDO
*
*   Read orbital information from Read Orbitals File; write summary
*   to  .dbg  file if option set
*
      IF (LDBPR(3) .AND. myid .EQ. 0) WRITE (99,300)
      NWIN = 0
    3 CONTINUE

      IF (myid .EQ. 0) READ (23,IOSTAT = IOS) NPY,NAKY,EY,MY
      CALL MPI_Bcast (IOS,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (NPY,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (NAKY,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (MY,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (EY,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      IF (IOS .EQ. 0) THEN
         CALL ALLOC (PNTRPA,MY,8)
         CALL ALLOC (PNTRQA,MY,8)
         CALL ALLOC (PNTRRA,MY,8)

         IF (myid .EQ. 0) THEN
            READ (23) PZY,(PA(I),I = 1,MY),(QA(I),I =1 ,MY)
            READ (23) (RA(I),I = 1,MY)
         ENDIF
         CALL MPI_Bcast (PZY,1,MPI_DOUBLE_PRECISION,0,
     &                         MPI_COMM_WORLD,ierr)
         CALL MPI_Bcast (PA,MY,MPI_DOUBLE_PRECISION,0,
     &                         MPI_COMM_WORLD,ierr)
         CALL MPI_Bcast (QA,MY,MPI_DOUBLE_PRECISION,0,
     &                         MPI_COMM_WORLD,ierr)
         CALL MPI_Bcast (RA,MY,MPI_DOUBLE_PRECISION,0,
     &                         MPI_COMM_WORLD,ierr)

         DO J = 1, NW
            IF ( (E(J) .LT. 0.D0) .AND.
     :           (NPY .EQ. NP(J)) .AND. (NAKY .EQ. NAK(J)) ) THEN
               PZ(J) = PZY
               E(J) = EY
               CALL INTRPQ (PA, QA, MY, RA, J, DNORM)
               IF (LDBPR(3) .AND. myid .EQ. 0) 
     &               WRITE (99,301) NP(J), NH(J), E(J), DNORM
               NWIN = NWIN + 1
            ENDIF
         ENDDO

         CALL DALLOC (PNTRPA)
         CALL DALLOC (PNTRQA)
         CALL DALLOC (PNTRRA)
         GOTO 3
      ENDIF
      IF (LDBPR(3) .AND. myid .EQ. 0)
     &      WRITE (99,*) ' orbitals renormalised;'
*
*   Return with an error code if all orbitals are not known
*
      IF (NWIN .LT. NW) THEN
         ierror = 1
         RETURN
      ENDIF
*
*   Schmidt orthogonalise the orbitals
*
      CALL ORTHSC

      IF (LDBPR(3) .AND. myid .EQ. 0)
     &      WRITE (99,*) 'orbitals orthogonalised and renormalised;'
      ierror = 0

  300 FORMAT (/'From SUBROUTINE LODRWF:'
     :        /' Orbital',8X,'Eigenvalue',19X,'Norm')
  301 FORMAT (2X,I2,A2,4X,1P,1D22.15,4X,1D22.15)

      RETURN
      END
