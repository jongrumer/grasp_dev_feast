************************************************************************
*                                                                      *
      SUBROUTINE MCP (RESTRT)
*                                                                      *
*   This routine controls the computation  and storage of the values   *
*   and all indices of the angular coefficients                        *
*                                                                      *
*                                       k                              *
*                   T  (ab)            V  (abcd)                       *
*                    rs                 rs                             *
*                                                                      *
*   k is the multipolarity of a two-particle Coulomb integral. a, b,   *
*   c and d are orbital sequence numbers.  r and s are configuration   *
*   state function indices.                                            *
*                                                                      *
*   Call(s) to: [LIB92]: ALCBUF, ALLOC, CONVRT, DALLOC, RKCO_GG,       *
*                        TNSRJJ.                                       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 21 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)

      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)

Cww      INTEGER PNJCUP,PNTJQS,PNTRIQ
      POINTER(PNTRIQ,RIQDUMMY(*))
      POINTER(PNJCUP,JCUPDUMMY(*))
      POINTER(PNTJQS,JQSDUMMY(*))

      LOGICAL DIAG,F0INT,LDBPA,LFORDR,LINCR,RESTRT
      CHARACTER*24 NAME
      CHARACTER*20 CNUM
      CHARACTER*2 CK,NH
      CHARACTER*500 CSF1
      CHARACTER*500 LINE1
*
      DIMENSION TSHELL(NNNW)
*
      POINTER (PLISTV,LLISTV(0:*)),(PNCFI,NCFI(*))
*
      POINTER (PLABEL,LABEL(6,*)),(PCOEFF,COEFF(*))
*
CGG      EXTERNAL COR,CORD
      EXTERNAL CORD
*
      COMMON/BUFFER/NBDIM,PLABEL,PCOEFF,NVCOEF
     :      /DEBUGA/LDBPA(5)
     :      /DEBUG/IBUG1,IBUG2,IBUG3,IBUG4,IBUG5,IBUG6
     :      /FOPARM/ICCUT
      COMMON/IDENT/NPOS(1000),NMR,NBEGIN2
     :      /MCPA/KMAX
     :      /MCPB/DIAG,LFORDR
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /STAT/PNTJQS,PNJCUP
*
*   Set the encoding key and its square
*
CGG begin
      PARAMETER (KEY = KEYORB)
CGG end
      PARAMETER (KEYSQ = KEY*KEY)
CGG      PARAMETER (KEY = 121, KEYSQ = 121*121)
*
*   Establish the cutoff criterion
*
      PARAMETER (CUTOFF = 1.0D-10)
*
*   Allocate storage to the array that stores the lengths
*   of the lists of V coefficients of each multipolarity
*
c     CALL ALLOC (PLISTV,KMAX+1,4)
      CALL ALLOC (PNCFI,NCF,4)
*
      DO K = 1,NCF
        NCFI(K) = 0
      ENDDO
*
*   Set the rank (zero) and parity (even) for the one-particle
*   coefficients
*
      KA = 0
      IOPAR = 1
*
*   INCOR is 0: only coefficients between open subshells required
*
      INCOR = 0
*
*   Allocate storage for the arrays in BUFFER
*
      CALL ALCBUF (1)
*      write(*,*) 'Hej2'
*
*   JA and JB respectively refer to the initial and final states
*   in the list of NCF configurations
*
      DO JJA = 1,NMR
        JA = NPOS(JJA)
        NCFI(JA) = 1
        DO JB = 1,NCF
*
          IF (JB .NE. JA) THEN
*
*   Call the MCT package to compute T coefficients
*
           CALL ONESCALAR(JA,JB,IA,IB,TSHELL)
CGG            CALL TNSRJJ (KA,IOPAR,JA,JB,IA,IB,TSHELL)
*
*   Write T coefficients that have magnitudes greater than the
*   cutoff criterion
*
            IF (IA .NE. 0) THEN
              IF (IA .NE. IB) THEN
                TCOEFF = TSHELL(1)
                IF (ABS (TCOEFF) .GT. CUTOFF) NCFI(JB) = 1
              ENDIF
            ENDIF
*
          ENDIF
*
*   Initialize
*
          NVCOEF = 0
*
*   Call the MCP package to generate V coefficients; ac and bd
*   are the density pairs
*
CGG          CALL RKCO (JA,JB,COR,CORD,INCOR)
          CALL RKCO_GG (JA,JB,CORD,0,1)
*
          DO I = 1,NVCOEF
            VCOEFF = COEFF(I)
            IF (ABS (VCOEFF) .GT. CUTOFF) THEN
              IA = LABEL(1,I)
              IB = LABEL(2,I)
              IC = LABEL(3,I)
              ID = LABEL(4,I)
              K  = LABEL(5,I)
              F0INT = (K .EQ. 0) .AND.
     :                (IA .EQ. IC) .AND.
     :                (IB .EQ. ID)
              IF (.NOT. F0INT) NCFI(JB) = 1
            ENDIF
          ENDDO
*
        ENDDO
*
      ENDDO
*
*   Deallocate storage that is no longer required
*
      CALL DALLOC (PNTRIQ)
      CALL DALLOC (PNTJQS)
      CALL DALLOC (PNJCUP)
      CALL ALCBUF (3)
c     CALL DALLOC (PLISTV)

      NINTER = 0
      DO K = 1,NCF
        IF (NCFI(K).EQ.1) NINTER = NINTER + 1
      ENDDO
     
      OPEN (14,FILE='rcsl.inp',FORM='FORMATTED',STATUS='OLD')
      OPEN (15,FILE='rcsl.out',FORM='FORMATTED',STATUS='UNKNOWN')

      REWIND (14)
      REWIND (15)

      DO I = 1,NBEGIN2
        READ(14,300) LINE1
        K = 500
    8   IF (LINE1(K:K) .EQ. ' ') THEN
          K = K-1
          IF (K .GT. 1) GOTO 8 
        ENDIF
        WRITE(15,300) LINE1(1:K)
      ENDDO
      DO I = 1,NCF
        DO J = 1,3
          READ(14,300) CSF1
          IF (NCFI(I).EQ.1) THEN
            K = 500
   11       IF (CSF1(K:K) .EQ. ' ') THEN
              K = K-1
              IF (K .GT. 1) GOTO 11
            ENDIF
            WRITE(15,300) CSF1(1:K) 
          ENDIF
        ENDDO
      ENDDO

      CLOSE (14)
      CLOSE (15)
*
      RETURN
*
  300 FORMAT (A)
      END
