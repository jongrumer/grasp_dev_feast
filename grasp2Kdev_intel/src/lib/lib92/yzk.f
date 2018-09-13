************************************************************************
*                                                                      *
      SUBROUTINE YZK (K,I,J)
*                                                                      *
*   This subroutine evaluates Hartree Y- and Z-functions:              *
*                                                                      *
*               (K)            (K)           (K)                       *
*              Y   (I,J;r) =  Z   (I,J;r) + W   (I,J;r)                *
*                                                                      *
*   where                                                              *
*                                                                      *
*    (K)                                                               *
*   Z   (I,J;r) =  I ( (s/r)   (P (s)*P (s) + Q (s)*Q (s)) ; 0 - r )   *
*                                I     J       I     J                 *
*                                                                      *
*   and                                                                *
*                                                                      *
*    (K)                    K+1                                        *
*   W   (I,J;r) =  I ( (r/s)   (P (s)*P (s) + Q (s)*Q (s)) ; r -       *
*                                I     J       I     J    INFINITY )   *
*                                                                      *
*   where  I ( G(r,s) ; range )  denotes the integral of G(r,s) over   *
*   range  in  s .  The Y-function is tabulated in  COMMON/TATB/  in   *
*   array  TB , the Z-function in array tA .                           *
*                                                                      *
*   Written by Farid A Parpia, at Oxford   Last updated: 06 Oct 1992   *
*   Modified by Anders Ynnerman, at Vanderbilt         : 03 Feb 1994   *
*   Modified by Jacek Bieron     at Vanderbilt         : 08 Feb 1994   *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
*
      parameter (klimit = 20)
      parameter (klimit1 = klimit +1)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
      LOGICAL KCALC(0:klimit)
      DIMENSION RHOP(NNN1),RTTK(NNN1,klimit),RTTKM(NNN1,klimit),
     :RTTK1(NNN1,klimit),
     :          RTTKM1(NNN1,klimit),RM(NNN1),WK(NNN1),TEMP(NNN1),
     :YK(NNN1),ZK(NNN1)
*
      POINTER (PNTRPF,PF(NNNP,*)),(PNTRQF,QF(NNNP,*))
*
      COMMON/CNC5/CNC5C(2:5,2:4)
     :      /DEF4/ACCY,NSCF,NSIC,NSOLV
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /NCC/C1,C2,C3,C4
     :      /ORB1/E(NNNW),GAMA(NNNW)
     :      /TATB/TA(NNN1),TB(NNN1),MTP
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
*
      DATA KCALC/klimit1*.FALSE./
      SAVE KCALC,RTTK,RTTKM,RTTK1,RTTKM1,RM
      EQUIVALENCE (TA(1),ZK(1)),(TB(1),YK(1))
*
*                         K
*   For  K > 0  compute  R   and store in  RTTK
*
      IF (.NOT.KCALC(K)) THEN
        KCALC(K) = .TRUE.
        if (k.gt.klimit) then
          print *, ' increase klimit in yzk.f '
          stop
        endif
        IF (K .GT. 0) THEN
           DO 1 II = 2,N
              RTTK(II,K) = R(II)**K
              RTTKM(II,K) = 1.D0/RTTK(II,K)
              RTTK1(II,K) = RTTK(II,K)*R(II)
              RTTKM1(II,K) = 1.D0/RTTK1(II,K)
              RM(II) = 1.D0/R(II)
    1      CONTINUE
        ELSE
           DO II = 2,N
              RM(II) = 1.D0/R(II)
           END DO
        ENDIF
      ENDIF
*
*   Determine maximum tabulation point as location beyond which
*   RHOP  (see comment statements below) would be zero; determine
*   other important locations
*
      MTP = MIN (MF(I),MF(J))
      MTPP1 = MTP+1
      MTPP3 = MTP+3
      MTPP4 = MTP+4
*
*   Compute RP(s)*(P (s)*P (s)+Q (s)*Q (s)) and store in RHOP
*                   I     J     I     J
*
      DO 2 II = 2,MTP
         RHOP(II) = RP(II)*(PF(II,I)*PF(II,J)+QF(II,I)*QF(II,J))
    2 CONTINUE
*
*   Fill array TEMP with R**K * RHOP
*
      TEMP(1) = 0.0D 00
      IF (K .EQ. 0) THEN
         DO 3 II = 2,MTP
            TEMP(II) = RHOP(II)
    3    CONTINUE
      ELSE
         DO 4 II = 2,MTP
            TEMP(II) = RTTK(II,K)*RHOP(II)
    4    CONTINUE
      ENDIF
*
*   Set an additional four points to zero
*
      DO 5 II = MTPP1,MTPP4
         TEMP(II) = 0.0D 00
    5 CONTINUE
*
*                                     K
*   Compute the first few values of  R  * ZK  using semi-open
*   Newton-Cotes formulae
*
      ZK(1) = 0.0D 00
      DO 7 II = 2,4
         SUM = 0.0D 00
         DO 6 KK = 2,5
            SUM = SUM+CNC5C(KK,II)*TEMP(KK)
    6    CONTINUE
         ZK(II) = SUM
    7 CONTINUE
*                         K
*   Compute remainder of R  * ZK: march out to MTP+3; use closed
*   Newton-Cotes formula
*
      DO 8 II = 5,MTPP3
         RTMP = C1*(TEMP(II-4)+TEMP(II  ))
     :         +C2*(TEMP(II-3)+TEMP(II-1))
     :         +C3* TEMP(II-2)
         ZK(II) = ZK(II-4) + RTMP
    8 CONTINUE
*
*                                       K   (K)
*   Determine the asymptotic value of  R * Z
*
*                   (0)
*   Correction to  z   : in the manner of  C Froese Fischer,
*   The Hartree-Fock Method for Atoms, John Wiley & Sons,
*   New York, 1977, p 235.
*
      IF (K .EQ. 0) THEN
*
         IF (I .EQ. J) THEN
            ZKLIM = 1.0D 00
         ELSE
            ZKLIM = 0.0D 00
         ENDIF
*
         DO 10 KK = MTPP3,MTP,-1
            DIF = ZK(KK)-ZKLIM
            IF (ABS (DIF) .GT. ACCY) THEN
               DO 9 II = KK,2,-4
                  ZK(II) = ZK(II)-DIF
    9          CONTINUE
            ENDIF
   10    CONTINUE
*
      ELSE
*
         ZKLIM = ZK(MTPP3)
*
      ENDIF
*
*   Tabulate  ZK  for entire internal grid
*
      IF (K .EQ. 0) THEN
*
         DO 11 II = MTPP4,N
            ZK(II) = ZKLIM
   11    CONTINUE
*
      ELSE
*
         DO 12 II = 2,MTPP3
            ZK(II) = ZK(II)*RTTKM(II,K)
   12    CONTINUE
*
         DO 13 II = MTPP4,N
            ZK(II) = ZKLIM*RTTKM(II,K)
   13    CONTINUE
*
      ENDIF
*
*   Start array WK / R**(K+1)
*
      NP4 = N+4
      DO 14 II = NP4,MTPP1,-1
         WK(II) = 0.0D 00
   14 CONTINUE
*
*   Fill array TEMP with RHOP / R**(K+1) ; set TEMP(1) = 0
*   to avoid 0/0 case
*
      TEMP(1) = 0.0D 00
      IF (K .EQ. 0) THEN
         DO 16 II = 2,MTP
            TEMP(II) = RHOP(II)*RM(II)
   16    CONTINUE
      ELSE
         DO 17 II = 2,MTP
            TEMP(II) = RHOP(II)*RTTKM1(II,K)
   17    CONTINUE
      ENDIF
*
*   Compute remainder of WK / R**(K+1): march in to the origin
*
      DO 18 II = MTP,2,-1
         WK(II) = WK(II+4)+C1*(TEMP(II  )+TEMP(II+4))
     :                    +C2*(TEMP(II+1)+TEMP(II+3))
     :                    +C3*(TEMP(II+2))
   18 CONTINUE
      WK(1) = 0.0D 00
*
*   Compute WK
*
      IF (K .EQ. 0) THEN
         DO 19 II = 2,MTP
            WK(II) = WK(II)*R(II)
   19    CONTINUE
      ELSE
         DO 20 II = 2,MTP
            WK(II) = WK(II)*RTTK1(II,K)
   20    CONTINUE
      ENDIF
*
*   Assemble solution
*
      YK(1) = 0.0D 00
      DO 21 II = 2,N
         YK(II) = ZK(II)+WK(II)
   21 CONTINUE
*
      RETURN
      END
