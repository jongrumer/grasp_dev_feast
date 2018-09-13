************************************************************************
*                                                                      *
      SUBROUTINE ZKF (K,I,J)
*                                                                      *
*   This subroutine evaluates Hartree Z-functionals:                   *
*                                                                      *
*              (k)                     k                               *
*             Z   [f(r);r] =  I ( (s/r)   f(s) ; 0 - r )               *
*                                                                      *
*   where  I ( g(r,s) ; range )  denotes the integral of g(r,s) over   *
*   range  in  s .    The Z-functional is tabulated in  COMMON/TATB/   *
*   in array  TB . The f-function is assumed tabulted in array  TA .   *
*                                                                      *
*   Written by Farid A Parpia, at Oxford   Last updated: 14 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
*
      DIMENSION RHOP(NNN1),RTTK(NNN1),TEMP(NNN1),ZK(NNN1)
*
      COMMON/CNC5/CNC5C(2:5,2:4)
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /NCC/C1,C2,C3,C4
     :      /TATB/TA(NNN1),TB(NNN1),MTP
*
      EQUIVALENCE (TB(1),ZK(1))
*
*                         k
*   For  k > 0  compute  r   and store in  RTTK
*
      IF (K .GT. 0) THEN
         DO 1 II = 2,N
            RTTK(II) = R(II)**K
    1    CONTINUE
      ENDIF
*
*   MTP is fed in through  COMMON/TATB/
*
      MTPP1 = MTP+1
      MTPP3 = MTP+3
      MTPP4 = MTP+4
*
*   Compute  RP(S)*F(S)  and store it in  RHOP
*
      DO 2 II = 2,MTP
         RHOP(II) = RP(II)*TA(II)
    2 CONTINUE
*
*   Fill array TEMP with r**k * RHOP
*
      TEMP(1) = 0.0D 00
      IF (K .EQ. 0) THEN
         DO 3 II = 2,MTP
            TEMP(II) = RHOP(II)
    3    CONTINUE
      ELSE
         DO 4 II = 2,MTP
            TEMP(II) = RTTK(II)*RHOP(II)
    4    CONTINUE
      ENDIF
*
*   Set an additional four points to zero
*
      DO 5 II = MTPP1,MTPP4
         TEMP(II) = 0.0D 00
    5 CONTINUE
*
*                                     k
*   Compute the first few values of  r  * ZK  using semi-open
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
*                         k
*   Compute remainder of r  * ZK: march out to MTP+3
*
      DO 8 II = 5,MTPP3
         ZK(II) = ZK(II-4)+C1*(TEMP(II-4)+TEMP(II  ))
     :                    +C2*(TEMP(II-3)+TEMP(II-1))
     :                    +C3* TEMP(II-2)
    8 CONTINUE
*                                       k   (k)
*   Determine the asymptotic value of  r * Z
*
*   Compute ZK
*
      ZKLIM = ZK(MTPP3)
*
      IF (K .EQ. 0) THEN
*
         DO 9 II = MTPP4,N
            ZK(II) = ZKLIM
    9    CONTINUE
*
      ELSE
*
         DO 10 II = 2,MTPP3
            ZK(II) = ZK(II)/RTTK(II)
   10    CONTINUE
*
         DO 11 II = MTPP4,N
            ZK(II) = ZKLIM/RTTK(II)
   11    CONTINUE
*
      ENDIF
*
      RETURN
      END
