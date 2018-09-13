************************************************************************
*                                                                      *
      SUBROUTINE GETINF
*                                                                      *
*   Interactively determines data useful for generating estimates of   *
*   the subshell radial wavefunctions.                                 *
*                                                                      *
*   Call(s) to: [LIB92]: GETYN, NUCPOT, RADGRD, SETQIC.                *
*               [RCI92]: SETISO.                                       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 21 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
      COMMON/iounit/istdi,istdo,istde

      LOGICAL GETYN,YES
*
      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /DEF2/C
     :      /DEF4/ACCY,NSCF,NSIC,NSOLV
     :      /DEF9/CVAC,PI
     :      /DEFAULT/NDEF
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /NPAR/PARM(2),NPARM
*
*   Open, check, load data from, and close the  .iso  file
*
      CALL SETISO ('isodata') ! NOW ALSO READS GRID PARAMETERS FROM ISODATA

*
*   Set defaults
*
      C = CVAC
*
!      IF (NPARM .EQ. 0) THEN ! Model nucleus as point source!
!         RNT = EXP (-65.0D 00/16.0D 00) / Z
!         H = 0.5D 00**4
!         N = MIN (220,NNNP)
!      ELSE              !
!         RNT = 2.0D-06  ! These are now read from isodata instead
!         H = 5.0D-02    !
!         N = NNNP       !
!      ENDIF
!      HP = 0.0D 00      !

*
*     Print grid parameters to screen
*     
      PRINT*
      PRINT*, 'DEFAULT GRID PARAMETERS'
      PRINT*, '----------------------------------------------------'
      PRINT*, 'RNT =', RNT 
      PRINT*, 'H   =', H 
      PRINT*, 'HP  =', HP
      PRINT*, 'N   =', N
      PRINT*, '----------------------------------------------------'
      PRINT*
*
      IF (NDEF.NE.0) THEN
         WRITE(istde,*) 'Change the default speed of light '
     &                , 'or radial grid parameters?'
*
         YES = GETYN ()
         IF (YES) THEN
*
*   Modify the speed of light
*
            WRITE(istde,*) 'The physical speed of light in '
     %             , 'atomic units is',CVAC, ';'
            WRITE(istde,*) 'revise this value?'
            YES = GETYN ()
            IF (YES) THEN
               WRITE(istde,*) 'Enter the revised value:'
               READ *,C
            ENDIF
*
*   Modify the parameters controlling the radial grid
*
            WRITE(istde,*) 'The default radial grid parameters '
     &                    , 'for this case are:'
            WRITE(istde,*) ' RNT = ', RNT
            WRITE(istde,*) ' H   = ', H
            WRITE(istde,*) ' HP  = ', HP
            WRITE(istde,*) ' N   = ', N
            WRITE(istde,*) 'revise these values?'
            YES = GETYN ()
            IF (YES) THEN
               WRITE(istde,*) 'Enter RNT:'
               READ *, RNT
               WRITE(istde,*) 'Enter H:'
               READ *, H
               WRITE(istde,*) 'Enter HP:'
               READ *, HP
               WRITE(istde,*) 'Enter N:'
               READ *, N
            ENDIF
*
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
      RETURN
      END
