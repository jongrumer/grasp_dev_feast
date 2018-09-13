************************************************************************
*                                                                      *
      FUNCTION FCO (K,IR,IA,IB)
*                                                                      *
*   This routine evaluates a coefficient                               *
*                                                                      *
*                                K                                     *
*                               f   (IA,IB)                            *
*                                IR                                    *
*                                                                      *
*   Here  K  is the multipolarity, IR  is the sequence number of the   *
*   configuration, and  IA  and  IB  are orbital  sequence  numbers.   *
*   ( I P Grant,  B J McKenzie,  P H Norrington, D F Mayers, and N C   *
*   Pyper,  Computer Phys Commun 21 (1980) 207-231, Eqs (6). )         *
*                                                                      *
*   Call(s) to: [LIB92]: CLRX, IQ.                                     *
*                                                                      *
*   Written by Farid A Parpia, at Oxford  Last revision: 21 Dec 1992   *
*                                                                      *
************************************************************************
! . Two calls to IQ are physically inlined here to reduce the overhead.
! . In the inlined text, the checks for the integer range (ref. 
!    lib92/iq.f for details) are removed since FCO is called only by
!    SETCOF and SETHAM where the input arguments always fall in the
!    suitable range.
!XHH 1997.03.05 
*
      IMPLICIT REAL*8          (A-H,O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
      LOGICAL LDBPA
      CHARACTER*2 NH
*
      COMMON/DEBUGA/LDBPA(5)
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
     :      /ORB10/NH(NNNW)

!XHH
CGG      PARAMETER (NNNWP = 30)
      INTEGER*4 iIQA
      POINTER (PNTRIQ,iIQA(NNNWP,1))
      COMMON/ORB2/NCF,NW,PNTRIQ
!
      IF (IA .EQ. IB) THEN
*
!         IQA = IQ (IA,IR)
         iqa = IBITS (iiqa((ia-1)/4+1,ir),8*MOD(ia-1,4),8)
*
         IF (K .EQ. 0) THEN
            FCO = DBLE ((IQA*(IQA-1))/2)
         ELSE
            IQF = NKJ(IA)+1
            IF (IQA .EQ. IQF) THEN
               KAPPA = NAK(IA)
               FAC = CLRX (KAPPA,K,KAPPA)*DBLE (IQA)
               FCO = -0.5D0*FAC*FAC
            ELSE
               FCO = 0.D0
            ENDIF
         ENDIF
*
      ELSE
*
         IF (K .EQ. 0) THEN
!            FCO = DBLE (IQ (IA,IR)*IQ (IB,IR))
            FCO = DBLE(IBITS(iiqa((ia-1)/4+1,ir),8*MOD(ia-1,4),8)
     &                *IBITS(iiqa((ib-1)/4+1,ir),8*MOD(ib-1,4),8))
         ELSE
            FCO = 0.D0
         ENDIF
*
      ENDIF
*
*   Debug printout
*
!      IF (LDBPA(3) .AND. (ABS (FCO) .GT. 0.0D 00))
!     :      WRITE (99,300) K,NP(IA),NH(IA),NP(IB),NH(IB),FCO,IR
!*
      RETURN
!*
!  300 FORMAT (/'  ',1I2
!     :        /' f    (',1I2,1A2,',',1I2,1A2,') = ',1PD21.14,
!     :        /'  ',1I3/)
!*
      END
