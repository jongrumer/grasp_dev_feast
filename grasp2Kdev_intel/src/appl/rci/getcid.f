************************************************************************
*                                                                      *
      SUBROUTINE GETCID (isofile, rwffile, idblk)
*                                                                      *
*   Interactively determines the data governing the CI problem.        *
*   iccut is replaced by an array iccutblk(1:nblock)
*                                                                      *
*   Call(s) to: [LIB92]: GETYN, NUCPOT, RADGRD, SETQIC.                *
*               [RCI92]: SETISO, SETRWF.                               *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 15 Dec 1992   *
*   Block version by Xinghong He          Last revision: 15 Jun 1998   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      CHARACTER*(*) isofile, rwffile, idblk(*)*8

      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)

      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL GETYN,LFORDR,LTRANS,LVP,LSE,LNMS,LSMS,YES,lshift
*
      POINTER (PNTRPF,PF(NNNP,*))
      POINTER (PNTRQF,QF(NNNP,*))
*
      ! set in setiso/lodiso
      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /NPAR/PARM(2),NPARM
     :      /NSMDAT/SQN,DMOMNM,QMOMB

      COMMON/DECIDE/LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
     :      /DEF2/C
     :      /DEF4/ACCY,NSCF,NSIC,NSOLV
     :      /DEF9/CVAC,PI
     :      /DEFAULT/NDEF
     :      /FOPARM/ICCUT
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /NPOT/ZZ(NNNP),NNUC
     :      /ORB1/E(NNNW),GAMA(NNNW)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
     :      /WFAC/WFACT
     :      /BLIM/IPRERUN,NCSFPRE,COEFFCUT1,COEFFCUT2
     :      /WHERE/IMCDF
     :      /QEDCUT/NQEDCUT,NQEDMAX

      POINTER (pncfblk, ncfblk(0:*))
      COMMON/hblock/nblock, pncfblk

      POINTER (pnevblk, nevblk(*))
      POINTER (pncmaxblk, ncmaxblk(*))
      COMMON/hblock2/pnevblk, pncmaxblk

      POINTER (piccutblk, iccutblk(*))
      COMMON/iccu/piccutblk

      COMMON/iounit/istdi,istdo,istde
      COMMON/hamshiftj/nshiftj(100),nasfshift(100,100),
     :       asfenergy(100,100),lshift


!
! Open, check, load data from, and close the  .iso  file
! Data loaded are: EMN,Z,/NPAR/,/NSMDAT/ where /.../ means whole 
! common block.
!
C      PRINT *, 'Calling SETISO ...'
      CALL SETISO (isofile)
!
! The speed of light, if non-default then spread from node-0
! Quantities to be obtained: C
!
      IF (NDEF .NE. 0) THEN
         WRITE (istde,*) 'Revise the physical speed of light (',CVAC,
     &                     ' in a.u.) ?'
         YES = GETYN ()
         IF (YES) THEN
            WRITE (istde,*) 'Enter the revised value:'
            READ *,C
         ELSE
            C = CVAC
         ENDIF
      ELSE
         C = CVAC
      ENDIF
!
! Treat some CSF's as perturbation ? Broadcasting is used only in 
! non-default mode, as the above case for speed of light.
! Quantities to be obtained: LFORDR
!
      IF (NDEF .NE. 0) THEN
         WRITE (istde,*) 'Treat contributions of some CSFs'
     &,                    ' as first-order perturbations?'
         LFORDR = GETYN ()
      ELSE
         LFORDR = .FALSE.
      ENDIF

! Get iccutblk() from the user-input

      IF (.NOT. LFORDR) THEN
         !...Default first
         DO i = 1, nblock
           iccutblk(i) = ncfblk(i)
         ENDDO
      ELSE

      ! Let master do the i/o, then broadcast
            WRITE (istde,*) 'There are ', nblock, 'blocks. They are:'
            WRITE (istde,*) '  block     J Parity     No of CSFs'
            DO i = 1, nblock
               WRITE (istde,*) i, idblk(i)(1:5), ncfblk(i)
            ENDDO

            WRITE (istde,*)
            WRITE (istde,*) 'Enter iccut for each block'
            DO jblock = 1, nblock
               WRITE (istde,*) 'Block ', jblock, '   ncf = ',
     &                        ncfblk(jblock)
     &                , ' id = ', idblk(jblock)(1:5)
  123          READ (istdi,*) ntmp
               IF (ntmp .GE. 0 .AND. ntmp .LE. ncfblk(jblock)) THEN
                  iccutblk(jblock) = ntmp
               ELSE
                  WRITE (istde,*) 'ICCUT out of range, re-enter:'
                  GOTO 123
               ENDIF
            ENDDO
      ENDIF

******************************************************************
!
! Pre-run ?
!
!     IF (IPRERUN .EQ. 0) THEN

!        WRITE (istde,*) ' Prerun with limited interaction?'
!        YES = GETYN ()

!        IF (YES) THEN
!           IPRERUN = 1
!           LTRANS = .FALSE.
!           LVP = .FALSE.
!           LNMS = .FALSE.
!           LSMS = .FALSE.
!           LSE = .FALSE.

!           WRITE (istde,*)  ' Give CSL cut'
!           READ *, NCSFPRE
!           WRITE (istde,*)  ' Give coefficient cut for H_0'
!           READ *, COEFFCUT1
!           WRITE (istde,*) ' Give coefficient cut for the transvers'
!    &,                  ' interaction'
!           GOTO 99
!        ENDIF
!     ENDIF
******************************************************************
!
! Include transverse ?
! Quantities to be obtained: LTRANS, WFACT
!
      WRITE (istde,*) 'Include contribution of H (Transverse)?'
      LTRANS = GETYN ()
      WRITE (istde,*) 'Modify all transverse photon frequencies?'
      YES = GETYN ()
      IF (YES) THEN
         WRITE (istde,*) 'Enter the scale factor:'
         READ *, WFACT
      ELSE
         WFACT = 1.0D 00
      ENDIF
!
! Other interactions ? One logical for each case. Done altogether
!
      WRITE (istde,*) 'Include H (Vacuum Polarisation)?'
      LVP = GETYN ()

      WRITE (istde,*) 'Include H (Normal Mass Shift)?'
      LNMS = GETYN ()

      WRITE (istde,*) 'Include H (Specific Mass Shift)?'
      LSMS = GETYN ()

      WRITE (istde,*) 'Estimate self-energy?'
      LSE = GETYN ()
      IF (LSE.EQV..TRUE.) THEN
         NQEDCUT = 1
         WRITE (istde,*)
     :'Largest n quantum number for including self-energy for orbital'
         WRITE (istde,*) 'n should be less or equal 8'
         READ *, NQEDMAX
      ELSE
         NQEDCUT = 0
      END IF

      IF (NDEF .EQ. 0) THEN
         lshift = .FALSE.
      ELSE
         WRITE (istde,*) 'Shift diagonal matrix elements?'
         lshift = GETYN ()
      END IF

      IF (NDEF.EQ.0) THEN
       IF (LTRANS) THEN
        WRITE(734,'(a)') 'y            ! Contribution of H Transverse?'
       ELSE
        WRITE(734,'(a)') 'n            ! Contribution of H Transverse?'
       END IF
       IF (YES) THEN
        WRITE(734,'(a)') 'y            ! Modify photon frequencies?'
        WRITE(734,*) WFACT,'! Scale factor'
       ELSE
        WRITE(734,'(a)') 'n            ! Modify photon frequencies?'
       END IF
       IF (LVP) THEN
        WRITE(734,'(a)') 'y            ! Vacuum polarization?'
       ELSE
        WRITE(734,'(a)') 'n            ! Vacuum polarization?'
       END IF
       IF (LNMS) THEN
        WRITE(734,'(a)') 'y            ! Normal mass shift?'
       ELSE
        WRITE(734,'(a)') 'n            ! Normal mass shift?'
       END IF
       IF (LSMS) THEN
        WRITE(734,'(a)') 'y            ! Specific mass shift?'
       ELSE
        WRITE(734,'(a)') 'n            ! Specific mass shift?'
       END IF
       IF (LSE) THEN
        WRITE(734,'(a)') 'y            ! Self energy?'
        WRITE(734,*) NQEDMAX, '! Max n for including self energy'
       ELSE
        WRITE(734,'(a)') 'n            ! Self energy?'
       END IF

      END IF

!
! Parameters controlling the radial grid
!
!  99  IF (NPARM .EQ. 0) THEN
!         RNT = EXP (-65.D0/16.D0) / Z
!         H = 0.5D0**4
!         N = MIN (220,NNNP)
!      ELSE            ! Already defined by setiso
!         RNT = 2.D-6  ! Jon Grumer (Lund, 2013)     
!         H = 5.D-2    !
!         N = NNNP     !
!      ENDIF
!      HP = 0.D0       !
*
*     Print grid parameters to screen
*     
C      PRINT*
C      PRINT*, 'DEFAULT GRID PARAMETERS'
C      PRINT*, '----------------------------------------------------'
C      PRINT*, 'RNT =', RNT 
C      PRINT*, 'H   =', H 
C      PRINT*, 'HP  =', HP
C      PRINT*, 'N   =', N
C      PRINT*, '----------------------------------------------------'
C      PRINT*
*      
      IF ( NDEF.NE.0) THEN
         WRITE (istde,*) 'The default radial grid parameters'
     &,                 ' for this case are:'
         WRITE (istde,*) ' RNT = ',RNT,';'
         WRITE (istde,*) ' H = ',H,';'
         WRITE (istde,*) ' HP = ',HP,';'
         WRITE (istde,*) ' N = ',N,';'
         WRITE (istde,*) ' revise these values?'
         YES = GETYN ()

         IF (YES) THEN
            WRITE (istde,*) 'Enter RNT:'
            READ *, RNT
            WRITE (istde,*) 'Enter H:'
            READ *, H
            WRITE (istde,*) 'Enter HP:'
            READ *, HP
            WRITE (istde,*) 'Enter N:'
            READ *, N
         ENDIF
      ENDIF
*
*   ACCY is an estimate of the accuracy of the numerical procedures
*
      ACCY = H**6
*
*   Set up the coefficients for the numerical procedures
*
      CALL SETQIC
*
*   Generate the radial grid and all associated arrays
*
      CALL RADGRD
*
*   Generate $- r \times V_ (r)$
*
      CALL NUCPOT
*
*   Load the radial wavefunctions
*
      CALL SETRWFA (rwffile)
*
*   Write the basic parameters of the model electron cloud to the
*   .res  file; this is the second record on the file --- the
*   first is the header (see SUBROUTINE SETRES)
*
      WRITE (imcdf) NELEC, NCF, NW, nblock ! ncf is ncftot, from setcsll
*
*   Write the nuclear parameters and $- r \times V_ (r)$
*   to the  .res  file; these are the third, fourth, and fifth
*   records on the file
*
      WRITE (imcdf) Z, EMN
      WRITE (imcdf) NPARM,(PARM(I), I = 1, NPARM)
      WRITE (imcdf) N, (ZZ(I), I = 1, N), NNUC
*
*   Write the physical effects specification to the  .res  file.
*   iccutblk() is now an array of length nblock.
*   This is the sixth record on the file
*
      WRITE (imcdf) C, LFORDR, (ICCUTblk(i), i = 1, nblock), 
     &              LTRANS, WFACT, LVP, LNMS, LSMS
*
*   Write the grid data to the  .res  file; this is the seventh
*   record on the file
*
      NP10 = N + 10
      WRITE (imcdf) RNT, H, HP, (R(I), I = 1, NP10),
     &           (RP(I), I = 1, NP10), (RPOR(I), I = 1, NP10) ! (imcdf)
*
*   Write out the interpolated radial wavefunctions; there are
*   2*NW such records; thus the total number of records written
*   at this point is 7+2*NW
*
      DO J = 1, NW
         WRITE (imcdf) E(J), GAMA(J), PZ(J), MF(J)
         WRITE (imcdf) (PF(I,J), I = 1, MF(J)), (QF(I,J), I = 1, MF(J))
      ENDDO
*
      RETURN
      END
