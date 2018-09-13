************************************************************************
*                                                                      *
      SUBROUTINE QED (jstate,SLFINT,UCF)
*                                                                      *
*   This  routine estimates corrections to the  energy levels due to   *
*   self-energy.                                                       *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, DALLOC, IQ, QUAD, RALLOC, SCREEN.      *
*               [RCI92]: FZALF, HOVLAP.                                *
*                                                                      *
*                                           Last update: 30 Oct 1992   *
*   Modified by Xinghong He                 Last update: 24 Jun 1997   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PNTRIQ,PNIVEC
      POINTER (PNTRIQ,RIQDUMMY)
      POINTER (PNIVEC,IVECDUMMY)
      CHARACTER*2 NH, npchar*1, nakchar
*
      DIMENSION UCF(*)
*
      POINTER (PNEVEC,EVEC(*))
      POINTER (PNTRPF,PF(NNNP,*))
      POINTER (PNTRQF,QF(NNNP,*))
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
!
! Pre-set tolerable number for iteration in finding effective
! nuclear charge. 
!
      maxiter = 20
*
*   Determine `generalised occupation numbers'
*
czou     CALL ALLOC (PNTUCF,NW,8)
!
! Modified so that UCFJ describes the current eigenstate
!
      DO 4 J = 1,NW
         UCFJ = 0.0D 00
!         DO 3 I = 1,NVEC
          I = jstate
            DO 2 II = 1,NCF
               UCFJ = UCFJ+DBLE (IQ (J,II))
     :                   *EVEC(II+(I-1)*NCF)**2
    2       CONTINUE
!    3    CONTINUE
c     print *, ucfj,'ucf'
         UCF(J) = UCFJ
c zou    UCF(J) = UCFJ/DBLE (NCF)
    4 CONTINUE
*
      DO 14 J = 1,NW
*
         NPJ = NP(J)
*
         IF (NPJ .LE. 8) THEN
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
*
            zeff=z
            ratio=  ratden(PTEMP,QTEMP,MFJ,NPJ,KAPPA,ZEFF)
            VALU =  ratio 
     :             *FZALF (NPJ,KAPPA,ZEFF)
     :             /DBLE (NPJ**3)
*
   13       SLFINT(J) = VALU*ZEFF**4/(PI*C**3)
c        print *, 'No. orb.=',j,' Zeff = ',zeff     
c    &   , 'Scale= ',ratio    
c    &   , 'S.E. = ',slfint(j)*2*13.6058,slfint(j)/ratio*2*13.6058
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
************************************************************************
*                                                                      *
      FUNCTION RATDEN (P,Q,MTPO,NP,KAPPA,Z)
*                                                                      *
*   This subprogram computes the overlap of the orbital tabulated in   *
*   the arrays  P  and  Q  with maximum tabulation point  MTPO  with   *
*   a hydrogenic orbital with parameters  NP  KAPPA  Z .               *
*                                                                      *
*   Call(s) to: [LIB92]: DCBSRW, QUAD.                                 *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last update: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
*
      DIMENSION P(NNNP),Q(NNNP)
*
      COMMON/GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /HORB/PH(NNNP),QH(NNNP)
     :      /TATB/TA(NNN1),TB(NNN1),MTP
*
*   Set up the hydrogenic orbital
*
      CALL DCBSRW (NP,KAPPA,Z,EH,PZH,PH,QH,MTPH)
*
*   Compute the overlap
*
      MTP = MIN (MTPH,MTPO)
      do i=2,mtp
       if (rp(i).le.0.0219) k=i
      enddo
      mtp=k
      TA(1) = 0.0D 00
      DO 1 I = 2,MTP
          TA(I) = (P(I)*P(I)+Q(I)*Q(I))*RP(I)
c         TA(I) = (P(I)*PH(I)+Q(I)*QH(I))*RP(I)
    1 CONTINUE
      CALL QUAD (RESULT)
      DO 2 I = 2,MTP
          TA(I) = (PH(I)*PH(I)+QH(I)*QH(I))*RP(I)
    2 CONTINUE
      CALL QUAD (RESULT1)
*
      RATDEN = RESULT/RESULT1
*
      RETURN
      END
