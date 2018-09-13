************************************************************************
*                                                                      *
      SUBROUTINE MCPIN (NAME,IK,NTESTG,INPCI)
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
*   Bug corrected 2005-10-18                                           *
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
      SAVE EVSC,NCFSC,NVECSC
      common/BLKTOT/NELECTOT,NCFTOT,NWTOT,NVECTOT,NVECSIZTOT,NBLOCKT
      common /BLK/NBLOCK,NCFBLK(10)



      NTESTL = 00
      NTEST = MAX(NTESTG,NTESTL)
      NTEST = 00

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

      REWIND (NF)
      DO 1000 IBLK = 1, NBLOCK
*
*   Read the CI vectors
* 
      CALL GETMIX(NAME,INPCI,IBLK)
*
*   Allocate memory for use in CITRAG
*
      CALL alloc (pscr, ncf*nvec, 8)
      CALL alloc (pciout, ncf*nvec, 8)
*
*   if data available on file read the datafile. Else sort the
*   MCP data into inegral based lists.
*
      READ(NF) NCFD,NWD,KAMAXD
      DO L = 1,KAMAX

*************
*
*. Offset for given L in shell matrices
*
          IF (L.EQ.1) THEN
            IIOFF = 1
          ELSE
            IIOFF = IIOFF + NSHL(L-1)** 2
          END IF
*  Corrected PER J

***************

          READ(NF) NINT,NCOEFF
          IF(NCOEFF*NINT.EQ.0) CYCLE
*
*   Allocate memory. Note that in this case we can,
*   since the data is sorted, allocate less 
*   memory than is done in qqsort.
*
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
**
**. Offset for given L in shell matrices
**
*        IF (L.EQ.1) THEN
*          IIOFF = 1
*        ELSE
*          IIOFF = IIOFF + NSHL(L-1)** 2
*        END IF
*   BUG THE OFFSET SHOULD BE BEFORE THE CYCLE STATEMENT!!
*
* Transform for given L
*
        IF (NSHL(L).GT.0) 
     &  CALL CITRAG(EVEC(1),NCF,NVEC,L,NSHL(L),
     &              CIROT(IIOFF),NINL(L),NTESTL,CIOUT, SCR )
*
*   Deallocate storage for this kappa.
*
        CALL DALLOC (PJANN)
        CALL DALLOC (PJBNN)
        CALL DALLOC (PPPINT)
        CALL DALLOC (PCNN)
        CALL DALLOC (PINTPTR)

      ENDDO

      CALL dalloc (pscr)
      CALL dalloc (pciout)
*
*   Write the rotated CI vectors on file
*
      IF(IBLK.EQ.1) THEN
        J = INDEX(NAME,' ')
        IF (INPCI.EQ.0) THEN
          OPEN (UNIT = 31,FILE=NAME(1:J-1)//'.cbm',FORM='UNFORMATTED',
     :    STATUS='UNKNOWN')
        ELSE
          OPEN (UNIT = 31,FILE=NAME(1:J-1)//'.bm',FORM='UNFORMATTED',
     :    STATUS='UNKNOWN')
        ENDIF

        WRITE(31) 'G92MIX'
        WRITE(31) NELEC,NCFTOT,NW,NVECTOT,NVECSIZE,NBLOCK
      ENDIF

      WRITE(31) IBLK,NCF,NVEC,IATJPO(1),IASPAR(1)
      WRITE(31) (IVEC(I),I = 1,NVEC)
c     WRITE(31) (IATJPO(I),IASPAR(I),I = 1,NVEC)
      WRITE(31) EAV,(EVAL(I),I = 1,NVEC)
      WRITE(31) ((EVEC(I+(J-1)*NCF),I = 1,NCF),J = 1,NVEC)

      CALL DALLOC (PNEVAL)
      CALL DALLOC (PNEVEC)
      CALL DALLOC (PNIVEC)
      CALL DALLOC (PIATJP)
      CALL DALLOC (PIASPA)
 1000 CONTINUE
      CLOSE (NF)
*
*   Close the  mixing  file
*
      CLOSE (30)
      CLOSE (31)
*
      RETURN
*
      END
