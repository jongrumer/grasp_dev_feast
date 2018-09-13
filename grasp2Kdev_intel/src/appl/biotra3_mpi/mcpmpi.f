************************************************************************
*                                                                      *
      SUBROUTINE MCP (NAME,startdir,IK,NTESTG,INPCI)
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

      PARAMETER (NLMAX = 40)
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
      CHARACTER*128 startdir
      CHARACTER*20 CNUM                                                 
      CHARACTER*2 CK,NH
*
      DIMENSION TSHELL(NNNW),LLISTT(NLMAX),NAKINV(NNNW),NSHL(NLMAX)
      DIMENSION NINL(NLMAX),SC(NVMAX)
*
C     NOT NEEDED
C      EXTERNAL COR,CORD,ITJPO
*
      POINTER (PNEVAL,EVAL(*))
      POINTER (PNEVEC,EVEC(*))
*P    POINTER (PNEVEC,EVEC(500))
      POINTER (PNIVEC,IVEC(*))
      POINTER (PIATJP,IATJPO(*)),(PIASPA,IASPAR(*))

      POINTER (PJANN,JANN(*)),(PPPINT,INTGRL(*)),
     :        (PJBNN,JBNN(*)),(PCNN,CNN(*)),(PINTPTR,INTPTR(*))

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
      COMMON /mpi/ myid, nprocs, ierr                   

* Locals
      POINTER (pscr, scr(*))
      POINTER (pciout, ciout(*))
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
      PARAMETER (KEYSQ = KEYORB*KEYORB)
CGG      PARAMETER (KEY = 121, KEYSQ = KEY*KEY)
*
*   Establish the cutoff criterion
*
      PARAMETER (CUTOFF = 1.0D-10)
*
      DATA KSTART/0/
      SAVE KSTART,EVSC,NCFSC,NVECSC

      KSTART = KSTART + 1


      NTESTL = 00
      NTEST = MAX(NTESTG,NTESTL)
      NTEST = 00
*
*   Read the CI vectors
* 
*R    print *, "mcp: NAME is ", trim(startdir)//'/'//trim(NAME)
      CALL GETMIX(trim(startdir)//'/'//NAME,INPCI)
*
*   Set up data for the initial and final state case respectively
*
      IF (IK.EQ.1) THEN
        DO I = 1,NLMAX
          NSHL(I) = NSHLII(I)
          NINL(I) = NINII(I)
        ENDDO

        DO I = 1,NNNW
          NAKINV(I) = NAKINVII(I)
        ENDDO

        DO J = 1,NLMAX
          DO I = 1,NLMAX
            NSHLP(I,J) = NSHLPII(I,J)
          ENDDO
        ENDDO

        DO J = 1,NNNW
          DO I = 1,NLMAX
            NSHLPP(I,J) = NSHLPPII(I,J)
          ENDDO
        ENDDO

        DO I = 1,20*NLMAX*NLMAX
          CIROT(I) = CICI(I)
        ENDDO

      ELSEIF (IK.EQ.2) THEN

        DO I = 1,NLMAX
          NSHL(I) = NSHLFF(I)
          NINL(I) = NINFF(I)
        ENDDO

        DO I = 1,NNNW
          NAKINV(I) = NAKINVFF(I)
        ENDDO

        DO J = 1,NLMAX
          DO I = 1,NLMAX
            NSHLP(I,J) = NSHLPFF(I,J)
          ENDDO
        ENDDO

        DO J = 1,NNNW
          DO I = 1,NLMAX
            NSHLPP(I,J) = NSHLPPFF(I,J)
          ENDDO
        ENDDO

        DO I = 1,20*NLMAX*NLMAX
          CIROT(I) = CFCI(I)
        ENDDO
      ENDIF
*
*   Check if angular data is available. If available read this data.
*   If not available calculate the data
*
      CALL ANGDATA(NAME,AVAIL,KAMAX)
      IF (.NOT.AVAIL) THEN
*
*   Open scratchfiles to dump the T coefficients for each kappa
*
        DO K = 1,KAMAX
Cww        OPEN (UNIT=40+K,STATUS='UNKNOWN',FORM='FORMATTED')
          OPEN (UNIT=80+K,STATUS='UNKNOWN',FORM='UNFORMATTED')
        ENDDO
*
*   Initialize the counters for the total number of T coefficients
*
        DO I = 1,NLMAX
          LLISTT(I) = 0
        ENDDO
*
*   Set the rank (zero) and parity (even) for the one-particle
*   coefficients
*
        KA = 0
        IOPAR = 1
*
*   JA and JB respectively refer to the initial and final states
*   in the list of NCF configurations
*
        DO JA = myid+1,NCF,nprocs
*R      DO JA = 1,NCF
          IF (MOD(JA,20).EQ.0.AND.KSTART.EQ.1) THEN 
            WRITE(*,*) ' JA1 =',JA
          ELSEIF (MOD(JA,20).EQ.0.AND.KSTART.EQ.2) THEN
            WRITE(*,*) ' JA2 =',JA
          ENDIF
*
          DO JB = 1,NCF
*
*   Call the MCT package to compute T coefficients
*
            IF (NTRANS.EQ.1) THEN
              COMP = .FALSE.
              IF (KSTART.EQ.1) THEN
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
Cww                  WRITE (40+NAKINV(IA),301) JA,JB,
Cww     :            NP(IA),NH(IA),NP(IA),NH(IA),TSHELL(IA)
                    WRITE (80+NAKINV(IA)) JA,JB,LAB,TSHELL(IA)
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
                    JAN = JA
                    JBN = JB 
                  ELSEIF (NORDII.EQ.1.AND.NORDFF.EQ.1) THEN
*
*   Experssion for reversed orbital ordering
*
                    LAB = IB*KEY + IA
                    JAN = JB
                    JBN = JA
                  ELSE                                                     
                    WRITE(*,*) 'SOMETHING WRONG'                           
                    STOP                                                   
                  ENDIF
Cww                WRITE (40+NAKINV(IA),301) JAN,JBN,
Cww     :          NP(IA),NH(IA),NP(IB),NH(IB),TSHELL(1)
                  WRITE (80+NAKINV(IA)) JAN,JBN,LAB,TSHELL(1)
*     :        JAN,JBN,LAB,TSHELL(1),'else write file',80+NAKINV(IA)
                ENDIF
              ENDIF
            ENDIF
            ENDIF
*
          ENDDO 
        ENDDO

      ENDIF

*
*   Deallocate storage that is no longer required. This was
*   allocated in lodcsl.
*
      CALL DALLOC (PNTRIQ)
      CALL DALLOC (PNTJQS)
      CALL DALLOC (PNJCUP)
*
*   Allocate memory for use in CITRAG
*
      CALL alloc (pscr, ncf*nvec, 8)
      CALL alloc (pciout, ncf*nvec, 8)
*
*   if data available on file read the datafile. Else sort the
*   MCP data into inegral based lists.
*
      DO L = 1,KAMAX
        IF (AVAIL) THEN
          READ(NF) NINT,NCOEFF
*
*   Allocate memory. Note that in this case we can,
*   since the data is sorted, allocate less 
*   memory than is done in qqsort.
*
          if (Ncoeff .GT. 0) then
              CALL ALLOC (PJANN,NCOEFF,4)                                       
              CALL ALLOC (PJBNN,NCOEFF,4)                                       
              CALL ALLOC (PPPINT,NINT,4)                                      
              CALL ALLOC (PCNN,NCOEFF,8)                                        
              CALL ALLOC (PINTPTR,NINT,4)
     
              DO I = 1,NINT                                                   
                READ(NF) INTGRL(I),INTPTR(I)                                 
              ENDDO                                                           
              DO I = 1,NCOEFF                                                 
                READ(NF) CNN(I),JANN(I),JBNN(I)                              
              ENDDO
          endif ! ncoeff > 0
        ELSE
*R        IF (LLISTT(L).GT.0) 
*R   :      CALL QQSORT (L,LLISTT(L),KSTART,NAME,KAMAX)
            CALL QQSORT (L,LLISTT(L),KSTART,NAME,KAMAX)

        ENDIF
*
*. Offset for given L in shell matrices
*
        IF (L.EQ.1) THEN
          IIOFF = 1
        ELSE
          IIOFF = IIOFF + NSHL(L-1)** 2
        END IF
*
* Transform for given L
*
        IF (NSHL(L).GT.0) 
     &  CALL CITRAG(EVEC(1),NCF,NVEC,L,NSHL(L),
     &              CIROT(IIOFF),NINL(L),NTESTL,CIOUT, SCR )
*
*   Deallocate storage for this kappa.
*
        if (Ncoeff .GT. 0) then
            CALL DALLOC (PJANN)
            CALL DALLOC (PJBNN)
            CALL DALLOC (PPPINT)
            CALL DALLOC (PCNN)
            CALL DALLOC (PINTPTR)
        endif !ncoeff > 0

      ENDDO


      CALL dalloc (pscr)
      CALL dalloc (pciout)
*
*   If angular data read from file close the file
*
      IF (AVAIL) CLOSE (NF)
*
*   Write the rotated CI vectors on file
*

      if (myid .eq. 0) then
          J = INDEX(NAME,' ')
          IF (INPCI.EQ.0) THEN
       OPEN (UNIT = 31,FILE=trim(startdir)//'/'//NAME(1:J-1)//'.cbm',
     : FORM='UNFORMATTED', STATUS='UNKNOWN')
          ELSE
       OPEN (UNIT = 31,FILE=trim(startdir)//'/'//NAME(1:J-1)//'.bm',
     : FORM='UNFORMATTED', STATUS='UNKNOWN')
          ENDIF
    
          WRITE(31) 'G92MIX'
          WRITE(31) NELEC,NCF,NW
          WRITE(31) NVEC
          WRITE(31) (IVEC(I),I = 1,NVEC)
          WRITE(31) (IATJPO(I),IASPAR(I),I = 1,NVEC)
          WRITE(31) EAV,(EVAL(I),I = 1,NVEC)
          WRITE(31) ((EVEC(I+(J-1)*NCF),I = 1,NCF),J = 1,NVEC)
*   Close the  mixing  file
*
          CLOSE (31)
       endif !myid=0


*Rasa start
*R     if(myid .eq. 0) then
*R     J = INDEX(NAME,' ')        
*R     OPEN (UNIT = 97,FILE=trim(startdir)//'/'//NAME(1:J-1)//'.fmt',
*R   : FORM='FORMATTED', STATUS='UNKNOWN')
*R     endif
*R    do j=1,NVEC
*R        print *, 'myid=',myid,' Rotated eigenvector: '
*R        do i=1,ncf
*R            WRITE(*,999) i, EVEC(I+(J-1)*NCF)
*R            if (myid .eq. 0)
*R                WRITE(97,999) i, EVEC(I+(J-1)*NCF)
*R        enddo
*R    enddo
*R    if (myid .eq. 0) close(97)
999   FORMAT(1X,I5,2X,1P,D19.10) 
*Rasa end
*      print '(4i10)', NELEC,NCF,NW,NVEC,(IVEC(I),I = 1,NVEC)
*     :      , (IATJPO(I),IASPAR(I),I = 1,NVEC)
*      print '(4f20.15)', ((EVEC(I+(J-1)*NCF),I = 1,NCF),J = 1,NVEC)
*
*
*   Program checking (assumes that the initial and final states are the same
*   but in different representation. If the CI rotation is correctly done then
*   the scalar product of the counter rotated CI coefficients should be 1 or -1
*
      IF (KSTART.EQ.1) THEN
        NCFSC = NCF
        NVECSC = NVEC
        IF (NCF*NVEC.LT.10000) THEN
          DO J = 1,NVEC
            DO I = 1,NCF
              EVSC(I+(J-1)*NCF) = EVEC(I+(J-1)*NCF)
            ENDDO
          ENDDO
        ENDIF
      ELSE
        IF (NCFSC.EQ.NCF.AND.NVECSC.EQ.NVEC.AND.NVEC.LE.NVMAX) THEN
          DO J = 1,NVEC
            SC(J) = 0.D0
            DO I = 1,NCF
              SC(J) = SC(J) + EVSC(I+(J-1)*NCF)*EVEC(I+(J-1)*NCF)
            ENDDO
          ENDDO
          WRITE(*,*) 
          WRITE(*,*) ' **************'
          WRITE(*,*) ' SCALAR PRODUCT'
          WRITE(*,*) ' **************'
          WRITE(*,*) 
          DO J = 1,NVEC
            WRITE(*,*) ' J, SC',J,SC(J)
          ENDDO
        ENDIF
      ENDIF
*
*   Deallocate storage for CI vectors. 
*   These arrays were allocated in getmix 
*
      CALL DALLOC (PNEVAL)
      CALL DALLOC (PNEVEC)
      CALL DALLOC (PNIVEC)
      CALL DALLOC (PIATJP)
      CALL DALLOC (PIASPA)
*
*  Close the angular files
*
      if (.NOT.AVAIL) then
         DO L = 1,KAMAX
Cww        CLOSE(L+40)
           CLOSE(L+80)
         ENDDO 
      endif
*
      RETURN
*
  300 FORMAT (/'From MCP:')
  301 FORMAT (' T_[',1I2,',',1I2,']',
     :   ' (',1I2,1A2,',',1I2,1A2,') = ',1PD19.12)
*
      END
