************************************************************************
*                                                                      *
      SUBROUTINE YPOT (J)
*                                                                      *
*   This subroutine  tabulates the  potential function Y(r) (Eq (14)   *
*   in  I P Grant, B J McKenzie, P H Norrington, D F Mayers, and N C   *
*   Pyper, Computer  Phys  Commun  21 (1980) 211) for orbital J. The   *
*   function is tabulated in the COMMON array YP.                      *
*                                                                      *
*   Call(s) to: [LIB92]: DRAW, YZK                                     *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last update: 17 Dec 1992   *
*   MPI version by Xinghong He            Last revision: 05 Aug 1998   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
      LOGICAL LDBPR
Cww      INTEGER PNTNDA,PNTNXA,PNTRDA,PNTRXA
      POINTER (PNTNDA,NDADUMMY), (PNTNXA,NXADUMMY)
      POINTER (PNTRDA,DADUMMY), (PNTRXA,XADUMMY)
      CHARACTER*2 NH
*
      POINTER (PNTRYA,YA(1))
      POINTER (PNTNYA,NYA(1))
*
      COMMON/DEBUGR/LDBPR(30)
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /NPOT/ZZ(NNNP),NNUC
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /POTE/YP(NNNP),XP(NNNP),XQ(NNNP)
     :      /SCF2/PNTRDA,PNTRXA,PNTRYA,
     :            PNTNDA,PNTNXA,PNTNYA,
     :            NDCOF,NXCOF,NYCOF,
     :            NDDIM,NXDIM,NYDIM
     :      /TATB/TA(NNN1),TB(NNN1),MTP
*
CGG      PARAMETER (KEYORB = 121)
      PARAMETER (KEY = KEYORB)

      !*** Need nprocs only
      COMMON /mpi/ myid, nprocs, ierr
!-----------------------------------------------------------------------
*   Debug printout: composition

      IF (LDBPR(29) .OR. LDBPR(30)) WRITE (99,300) NP(J),NH(J)
*
*   Initialize array YP with the nuclear potential piece
*
*   Since YA() below  contains contributions from THIS node only,
*   the initialization should be in consistence with that.

      DO I = 1, N
         YP(I) = ZZ(I) / nprocs
      ENDDO

      DO INDEX = 1, NYCOF

         ! Decode information in label
         LABEL = NYA(INDEX)
         K     = MOD (LABEL, KEY)
         LABEL = LABEL / KEY
         IOY1  = MOD (LABEL, KEY)
         IOY2  = LABEL / KEY
         COEFF = YA(INDEX)

         IF (LDBPR(29)) THEN
            WRITE (99,301) K,COEFF,NP(IOY1),NH(IOY1),NP(IOY2),NH(IOY2)
         ENDIF

         CALL YZK (K, IOY1, IOY2)     ! Accumulate contributions
         DO I = 1, N
            YP(I) = YP(I) - COEFF * TB(I)
         ENDDO
      ENDDO
*
*   Debug printout
*
      IF (LDBPR(30)) THEN
         WRITE (99,302)
         NB3 = N/3
         IF (3*NB3 .EQ. N) THEN
            NROWS = NB3
         ELSE
            NROWS = NB3+1
         ENDIF
         DO 4 II = 1,NROWS
            II1 = II
            II2 = II1+NROWS
            II3 = II2+NROWS
            IF (II3 .LE. N) THEN
               WRITE (99,303) R(II1),YP(II1),R(II2),YP(II2),
     :                        R(II3),YP(II3)
            ELSEIF (II2 .LE. N) THEN
               WRITE (99,303) R(II1),YP(II1),R(II2),YP(II2)
            ELSE
               WRITE (99,303) R(II1),YP(II1)
            ENDIF
    4    CONTINUE
         CALL DRAW (YP,1.0D 00,YP,0.0D 00,N)
      ENDIF
*
      RETURN
*
  300 FORMAT (//' Direct potential for ',1I2,1A2,' orbital :'//)
  301 FORMAT (/25X,'(',1I2,')'
     :        /1X,1PD21.14,'* Y    (',1I2,1A2,',',1I2,1A2,')')
  302 FORMAT (//3(' --------- r --------- ------- Y (r) -------'))
  303 FORMAT (1P,6(1X,1D21.14))
*
      END
