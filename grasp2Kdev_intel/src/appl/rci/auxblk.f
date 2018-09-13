      SUBROUTINE auxblk (j2max, atwinv)
      IMPLICIT REAL*8          (A-H, O-Z)

      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)

*      ...Need LTRANS,LVP,LNMS,LSMS
      LOGICAL LFORDR,LTRANS,LVP,LSE,LNMS,LSMS,LDBPG
      COMMON/DECIDE/LFORDR,LTRANS,LVP,LSE,LNMS,LSMS

*      ...Need EMN, N,RP, ZDIST, TB
      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     &      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /NCDIST/ZDIST(NNNP)
     :      /TATB/TA(NNN1),TB(NNN1),MTP

*      ...Initialize NCOEI,FRSTCO
      LOGICAL FRSTCO
      POINTER (PCOEIL,COEILDUMMY)
      POINTER (PCOEVL,COEVLDUMMY)
      COMMON/COEILS/NDCOEA,NCOEI,PCOEIL,PCOEVL,FRSTCO

*      ...Initialize NVPI,FRSTVP
      LOGICAL FRSTVP
      POINTER (PINDVP,INDVPDUMMY)
      POINTER (PVALVP,VALVPDUMMY)
      COMMON/VPILST/NDVPA,NVPI,PINDVP,PVALVP,FRSTVP

*      ...Initialize NTPI,FIRST which are to be used in setham-brintX
      LOGICAL FIRST
      POINTER (PINDT1,INDT1DUMMY)
      POINTER (PINDT2,INDT2DUMMY)
      POINTER (PINDT3,INDT3DUMMY)
      POINTER (PINDT4,INDT4DUMMY)
      POINTER (PINDT5,INDT5DUMMY)
      POINTER (PINDT6,INDT6DUMMY)
      POINTER (PVALT1,VALT1DUMMY)
      POINTER (PVALT2,VALT2DUMMY)
      POINTER (PVALT3,VALT3DUMMY)
      POINTER (PVALT4,VALT4DUMMY)
      POINTER (PVALT5,VALT5DUMMY)
      POINTER (PVALT6,VALT6DUMMY)
      COMMON/BILST/PINDT1,PINDT2,PINDT3,PINDT4,PINDT5,PINDT6,
     :             PVALT1,PVALT2,PVALT3,PVALT4,PVALT5,PVALT6,
     :             NDTPA(6),NTPI(6),FIRST(6)

*      ...Initialize NKEI,FRSTKI
      LOGICAL FRSTKI
      POINTER (PINDKE,INDKEDUMMY)
      POINTER (PVALKE,VALKEDUMMY)
      COMMON/KEILST/NDKEA,NKEI,PINDKE,PVALKE,FRSTKI

*      ...Initialize NVINTI,FRSTVI
      LOGICAL FRSTVI
      POINTER (PVINIL,VINILDUMMY)
      POINTER (PVINVL,VINVLDUMMY)
      COMMON/VINLST/NDVIN,NVINTI,PVINIL,PVINVL,FRSTVI

*      ...Need NW
      POINTER (PNTRIQ,RIQDUMMY)
      COMMON/ORB2/NCF,NW,PNTRIQ

      COMMON/iounit/istdi,istdo,istde

************************************************************************

      FRSTCO = .TRUE.
      NCOEI = 0

      IF (LTRANS) THEN
*        ...Check the maximum numbers of orbtitals allowed in brint.f
         i = NNNW
         IF (J2MAX.EQ.11) THEN
            i = 114
         ELSEIF (J2MAX .EQ. 12) THEN
            i = 112
         ELSEIF (J2MAX .EQ. 13) THEN
            i = 110
         ELSEIF (J2MAX .EQ. 14) THEN
            i = 108
         ELSEIF (J2MAX .EQ. 15) THEN
            i = 106
         ELSEIF (J2MAX .EQ. 16) THEN
            i = 105
         ELSEIF (J2MAX .EQ. 17) THEN
            i = 103
         ELSEIF (J2MAX .EQ. 18) THEN
            i = 101
         ELSEIF (J2MAX .EQ. 19) THEN
            i = 100
         ELSEIF (J2MAX .GT. 20) THEN
            i =  90
         ENDIF

         IF (i .LT. NW) THEN
            WRITE (istde,*) 'In setham. The number of orbitals is too'
            WRITE (istde,*) 'large for the brint routine'
            STOP
         ENDIF

         DO I = 1,6
            FIRST(I) = .TRUE.
            NTPI(I) = 0
         ENDDO
      ENDIF
*
*   Initialisations for the vacuum polarisation corrections
*
      IF (LVP) THEN
         CALL NCHARG
         CALL VACPOL
         DO K = 2,N
            ZDIST(K) = TB(K)*RP(K)
         ENDDO
         FRSTVP = .TRUE.
         NVPI = 0
      ENDIF
*
*   Initialisations for nuclear translational energy corrections
*
      IF (EMN .GT. 0.d0) THEN
         ATWINV = 1.D0/EMN
         IF (LNMS) THEN
            FRSTKI = .TRUE.
            NKEI = 0
         ENDIF
         IF (LSMS) THEN
            FRSTVI = .TRUE.
            NVINTI = 0
         ENDIF
      ELSE
         ! atwinv will not be used
         LNMS = .FALSE.
         LSMS = .FALSE.
      ENDIF

      RETURN
      END
