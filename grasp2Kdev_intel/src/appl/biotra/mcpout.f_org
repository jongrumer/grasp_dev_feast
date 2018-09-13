************************************************************************
*                                                                      *
      SUBROUTINE MCPOUT (NAME,IK,NTESTG,INPCI)
*                                                                      *
*   This routine controls the computation  and storage of the values   *
*   and all indices of the angular coefficients                        *
*                                                                      *
*                                                                      *
*                   T  (ab)                                            *
*                    rs                                                *
*                                                                      *
*   k is the multipolarity of a two-particle Coulomb integral. a, b,   *
*   are orbital sequence numbers.  r and s are configuration           *
*   state function indices.                                            *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, DALLOC, TNSRJJ                         *
*               [GENMCP]: QSORT.                                       *
*                                                                      *
*   Written by Per Jonsson                Last revision: JUne 1996     *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)

      PARAMETER (NLMAX = 20)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)                                            
CGG begin
      PARAMETER (KEY = KEYORB)
CGG end
      PARAMETER (NF = 200)
      PARAMETER (NVMAX = 100)

Cww      INTEGER PNJCUP,PNTJQS,PNTRIQ                                   
      POINTER (PNJCUP,JCUPDUMMY)                                        
      POINTER (PNTJQS,JQSDUMMY)                                         
      POINTER (PNTRIQ,RIQDUMMY)                                         

      LOGICAL DIAG,F0INT,LDBPA,LFORDR,LINCR,RESTRT,COMP,AVAIL
      CHARACTER*24 NAME
      CHARACTER*20 CNUM                                                 
      CHARACTER*2 CK,NH
*
      DIMENSION TSHELL(NNNW),LLISTT(NLMAX),NAKINV(NNNW),NSHL(NLMAX)
      DIMENSION NINL(NLMAX),SC(NVMAX)
*
      EXTERNAL COR,CORD,ITJPO
*
      POINTER (PNEVAL,EVAL(1))
      POINTER (PNEVEC,EVEC(500))
      POINTER (PNIVEC,IVEC(1))
      POINTER (PIATJP,IATJPO(1)),(PIASPA,IASPAR(1))

      POINTER (PJANN,JANN(1)),(PPPINT,INTGRL(1)),
     :        (PJBNN,JBNN(1)),(PCNN,CNN(1)),(PINTPTR,INTPTR(1))

      DIMENSION CIROT(20*NLMAX*NLMAX),EVSC(10000)

      COMMON/DEBUGA/LDBPA(5)
     :      /DEBUG/IBUG1,IBUG2,IBUG3,IBUG4,IBUG5,IBUG6
     :      /DEF1/EMN,IONCTY,NELEC,Z
     :      /EIGVAL/EAV,PNEVAL
     :      /EIGVEC/PNEVEC
     :      /FOPARM/ICCUT
     :      /MCPA/KMAX
     :      /MCPB/DIAG,LFORDR
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /PRNT/NVEC,PNIVEC,NVECMX
     :      /SYMA/PIATJP,PIASPA
     :      /STAT/PNTJQS,PNJCUP
     :      /JQJC/JQJ1,JQJ2,ITJQJ1(40),ITJQJ2(40),NTRANS
      COMMON/ORBORD/NORDII,NORDFF

      COMMON/MCPDATA/PJANN,PJBNN,PPPINT,PCNN,PINTPTR,NCOEFF,NINT

      COMMON/SBDAT/NAKINVII(NNNW),NSHLII(NLMAX),NSHLPII(NLMAX,NLMAX),   
     :             NAKINVFF(NNNW),NSHLFF(NLMAX),NSHLPFF(NLMAX,NLMAX),   
     :             NSHLPPII(NLMAX,NNNW),NSHLPPFF(NLMAX,NNNW),
     :             NINII(NLMAX),NINFF(NLMAX),IKAPPA(NLMAX),KAMAX

      COMMON/SBDAT1/NSHLP(NLMAX,NLMAX),NSHLPP(NLMAX,NNNW)

* Locals
      POINTER (pscr, scr(1))
      POINTER (pciout, ciout(1))
*
*  Initial state commons
*
      COMMON/CIIMAT/CICI(20*NLMAX*NLMAX)
*
*  Final state commons
*
      COMMON/CFFMAT/CFCI(20*NLMAX*NLMAX)

*
*   Set the encoding key and its square
*
      PARAMETER (KEYSQ = KEY*KEY)
CGG      PARAMETER (KEY = 121, KEYSQ = KEY*KEY)
*
*   Establish the cutoff criterion
*
      PARAMETER (CUTOFF = 1.0D-10)
*
      common /BLK/NBLOCK,NCFBLK(10)

      print *,'NBLOCK,(NCFBLK(i),i=1,NBLOCK)'
      print *,NBLOCK,(NCFBLK(i),i=1,NBLOCK)

      NTESTL = 00
      NTEST = MAX(NTESTG,NTESTL)
      NTEST = 00
*
*   Set the rank (zero) and parity (even) for the one-particle
*   coefficients
*
        KA = 0
        IOPAR = 1
*
*   Check if angular data is available. If available read this data.
*   If not available calculate the data
*
      CALL ANGDATA(NAME,AVAIL,KAMAX)
      print *,'AVAIL=',AVAIL
      IF (AVAIL) RETURN

      IF (IK.EQ.1) THEN
        DO I = 1,NNNW
          NAKINV(I) = NAKINVII(I)
        ENDDO
      ELSE
        DO I = 1,NNNW
          NAKINV(I) = NAKINVFF(I)
        ENDDO
      ENDIF



         print *, ' open sorted ang. file .TB(NF)', NF
          J = INDEX(NAME,' ')
          OPEN (UNIT = NF,FILE=NAME(1:J-1)//'.TB',
     :          STATUS='UNKNOWN',FORM='UNFORMATTED')

        NCF0 = 1
      DO 1000 IBLK = 1, NBLOCK
*
*   Open scratchfiles to dump the T coefficients for each kappa
*
        DO K = 1,KAMAX
          OPEN (UNIT=80+K,STATUS='UNKNOWN',FORM='UNFORMATTED')
        ENDDO
*
*   Initialize the counters for the total number of T coefficients
*
        DO I = 1,NLMAX
          LLISTT(I) = 0
        ENDDO
*
*   JA and JB respectively refer to the initial and final states
*   in the list of NCF configurations
*
        DO JA = NCF0,NCFBLK(IBLK)
          IF (MOD(JA,20).EQ.0.AND.IK.EQ.1) THEN 
            WRITE(*,*) ' JA1 =',JA,JA-NCF0+1
          ELSEIF (MOD(JA,20).EQ.0.AND.IK.EQ.2) THEN
            WRITE(*,*) ' JA2 =',JA,JA-NCF0+1
          ENDIF
*
          DO JB = NCF0,NCFBLK(IBLK)
*
*   Call the MCT package to compute T coefficients
*
            IF (NTRANS.EQ.1) THEN
              COMP = .FALSE.
              IF (IK.EQ.1) THEN
                DO IBB = 1,JQJ1
                  IF (ITJPO(JA).EQ.ITJQJ1(IBB)) COMP = .TRUE.
                ENDDO
              ELSE
                DO IBB = 1,JQJ2
                  IF (ITJPO(JA).EQ.ITJQJ2(IBB)) COMP = .TRUE.
                ENDDO
              ENDIF
            ELSE
              COMP = .TRUE.
            ENDIF


            IF (COMP) THEN
*            write(*,*) JA,JB
            CALL TNSRJJ (KA,IOPAR,JA,JB,IA,IB,TSHELL)
            IF (IA .NE. 0) THEN
              IF (IA .EQ. IB) THEN
                DO IA = 1,NW
*
*   If T coefficient is greater than zero and the kappa quantum numbers
*   of the two orbitals are the same dump to file
*   In a later version use a buffer with a reasonable record length
*
                  IF (DABS(TSHELL(IA)).GT.CUTOFF) THEN
                    LLISTT(NAKINV(IA)) = LLISTT(NAKINV(IA)) + 1
                       LAB = IA*KEY + IA
                    WRITE (80+NAKINV(IA)) JA-NCF0+1,JB-NCF0+1,
     &                     LAB,TSHELL(IA)
                  ENDIF
                ENDDO
              ELSE
                IF (DABS(TSHELL(1)).GT.CUTOFF.AND.
     :            NAK(IA).EQ.NAK(IB)) THEN
                  LLISTT(NAKINV(IA)) = LLISTT(NAKINV(IA)) + 1
                  IF (NORDII.EQ.0.AND.NORDFF.EQ.0) THEN 
*  
*   Experssion for normal orbital ordering
*
                    LAB = IA*KEY + IB                      
                    JAN = JA-NCF0+1
                    JBN = JB-NCF0+1 
                  ELSEIF (NORDII.EQ.1.AND.NORDFF.EQ.1) THEN
*
*   Experssion for reversed orbital ordering
*
                    LAB = IB*KEY + IA
                    JAN = JB-NCF0+1
                    JBN = JA-NCF0+1
                  ELSE                                                     
                    WRITE(*,*) 'SOMETHING WRONG'                           
                    STOP                                                   
                  ENDIF
                  WRITE (80+NAKINV(IA)) JAN,JBN,LAB,TSHELL(1)
                ENDIF
              ENDIF
            ENDIF
            ENDIF
*
          ENDDO 
        ENDDO


*
*   sort the MCP data into inegral based lists.
*
      DO L = 1,KAMAX
          IF (LLISTT(L).GT.0) THEN
            CALL QQSORT (L,LLISTT(L),IK,NAME,KAMAX)
          ELSE
            IF (L.EQ.1) WRITE(NF) NCF,NW,KAMAX
            WRITE(NF) LLISTT(L),LLISTT(L)
          ENDIF
      ENDDO
*
*  Close the angular files
*
      DO L = 1,KAMAX
        CLOSE(L+80,STATUS='DELETE')
      ENDDO 
*
      NCF0 =NCFBLK(IBLK)+1
 1000 CONTINUE
*
*   Deallocate storage that is no longer required. This was
*   allocated in lodcsl.
*
      CALL DALLOC (PNTRIQ)
      CALL DALLOC (PNTJQS)
      CALL DALLOC (PNJCUP)
      RETURN
      END
