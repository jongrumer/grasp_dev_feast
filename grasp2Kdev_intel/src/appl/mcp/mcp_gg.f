************************************************************************
*                                                                      *
      SUBROUTINE mcp (nb, RESTRT, myid, nprocs, fhead)
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
*               [GENMCP]: FNDBEG, SETSDA, SORT.                        *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 28 Sep 1993   *
*   Modified by C. Froese Fischer for block computation.               *
*   Modified by Gediminas Gaigalas for new spin-angular integration.   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
      POINTER (PNJCUP,JCUPDUMMY)
      POINTER (PNTJQS,JQSDUMMY)
      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL DIAG,F0INT,LDBPA,LFORDR,LINCR,RESTRT
      CHARACTER*2 NH, fhead*(*)

      DIMENSION TSHELL(NNNW)

      POINTER (PLISTV,LLISTV(0:*))

      POINTER (PLABEL,LABEL(6,*))
      POINTER (PCOEFF,COEFF(*))

CGG      EXTERNAL COR, CORD, outsda
      EXTERNAL CORD, outsda

      COMMON/BUFFER/NBDIM,PLABEL,PCOEFF,NVCOEF
     :      /DEBUGA/LDBPA(5)
     :      /DEBUG/IBUG1,IBUG2,IBUG3,IBUG4,IBUG5,IBUG6
     :      /DEFAULT/NDEF
     :      /FOPARM/ICCUT(100)
     :      /MCPA/KMAX
     :      /MCPB/DIAG,LFORDR
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /STAT/PNTJQS,PNJCUP

CGG      PARAMETER (KEYORB =121)                ! encoding key
      PARAMETER (KEY = KEYORB, KEYSQ = KEYORB*KEYORB)

! In order for serial and parallel programs to give the same result

       PARAMETER (CUTOFF = 1.0D-10)          ! cutoff criterion
c     PARAMETER (CUTOFF = 0.D0)              ! cutoff criterion
!-----------------------------------------------------------------------
!xhh - swap not used
!      OPEN (20,FILE='sms.20',FORM='UNFORMATTED',
!     :      STATUS='UNKNOWN')
*
*   Allocate storage to the array that stores the lengths
*   of the lists of V coefficients of each multipolarity
*
      CALL ALLOC (PLISTV, KMAX + 1, 4)
*
*   If this is a restart, determine the pair of CSFs with which
*   the computation should begin.
*   If not restart, initialize the counters for the total number of
*   T coefficients and V coefficients of each multipolarity.
*
      IF (RESTRT) THEN
         PRINT *, 'Not ready for RESTART'
         STOP
         CALL FNDBEG (JASTRT, JBSTRT, npos, LLISTT, LLISTV)
      ELSE
         JASTRT = 1 + myid
         npos = 0
         LLISTT = 0
         DO K = 0, KMAX
            LLISTV(K) = 0
         ENDDO
      ENDIF
*
*   Allocate storage for the arrays in BUFFER
*
      CALL ALCBUF (1)
*
*   Write header to  .dbg  file if appropriate
*
      IF (LDBPA(2) .OR. LDBPA(3)) WRITE (99,300)
*
*   JA and JB respectively refer to the initial and final states
*   in the list of NCF configurations
*
      DO 5 JA = JASTRT, NCF, nprocs
*
!xhh - Lower triangle by rows
!         IF (DIAG .OR. (LFORDR .AND. (JA .GT. ICCUT(nb))) ) THEN
!            JBSTRT = JA
!         ELSE
!            JBSTRT = 1
!         ENDIF

*
         JBSTRT = 1
         DO 4 JB = JBSTRT, JA


         IF (DIAG .OR. (LFORDR .AND. (JB .GT. ICCUT(nb))) ) THEN
            IF (JB.NE.JA) CYCLE
         END IF


*   LINCR is .TRUE. if npos is to be incremented by 1; there
*   is always a diagonal element in each column

            IF (JB .NE. JA) THEN
               LINCR = .TRUE.
            ELSE
               npos = npos + 1
               LINCR = .FALSE.
            ENDIF

            IF (JB .NE. JA) THEN
               ! Compute T coefficients
           CALL ONESCALAR(JA,JB,IA,IB,TSHELL)
CGG               CALL TNSRJJ (0, 1, JA, JB, IA, IB, TSHELL)

               ! Store it (TSHELL) if it's greater than CUTOFF
               IF (IA .NE. 0 .AND. IA .NE. IB .AND.
     &                 ABS (TSHELL(1)) .GT. CUTOFF) THEN
!xhh -  Debug output disabled
!                  IF (LDBPA(2)) WRITE (99,301) JB,JA,
!     :                  NP(IA),NH(IA),NP(IB),NH(IB),TSHELL(1)
                  npos = npos + 1
                  LINCR = .FALSE.
                  LLISTT = LLISTT + 1
                  LAB = MIN (IA,IB) * KEY + MAX (IA,IB)
                  WRITE (31) JA, npos, LAB, TSHELL(1)
               ENDIF
            ENDIF
*
*   Call the MCP package to generate V coefficients; ac and bd
*   are the density pairs. NVCOEF is initialized here but changed
*   in /rkco/cor[d] via COMMON.
*
            NVCOEF = 0
CGG            CALL RKCO (JA, JB, COR, CORD, 0)
            CALL RKCO_GG (JA, JB, CORD, 0, 1)
            DO 2 I = 1, NVCOEF
               VCOEFF = COEFF(I)
               IF (ABS (VCOEFF) .GT. CUTOFF) THEN
                  IA = LABEL(1,I)
                  IB = LABEL(2,I)
                  IC = LABEL(3,I)
                  ID = LABEL(4,I)
                  K  = LABEL(5,I)
                  F0INT = (K .EQ. 0) .AND.
     :                    (IA .EQ. IC) .AND.
     :                    (IB .EQ. ID)
                  IF (.NOT. F0INT) THEN
!xhh -  Debug output disabled
!                      IF (LDBPA(3)) WRITE (99,302) K,JB,JA,
!     :                  NP(IA),NH(IA),NP(IB),NH(IB),
!     :                  NP(IC),NH(IC),NP(ID),NH(ID),VCOEFF

* Swap index to make sure IA <= IC, IB <= ID and record the number
* of swaps
                     NSWAP = 0
                     IF (IA .GT. IC) THEN
                        ISWAP = IC
                        IC = IA
                        IA = ISWAP
                        NSWAP = NSWAP + 1
                     ENDIF
                     IF (IB .GT. ID) THEN
                        ISWAP = ID
                        ID = IB
                        IB = ISWAP
                        NSWAP = NSWAP + 1
                     ENDIF

                     IF (LINCR) THEN
                        npos = npos + 1
                        LINCR = .FALSE.
                     ENDIF

                     LLISTV(K) = LLISTV(K) + 1
                     LAC = IA * KEY + IC
                     LBD = IB * KEY + ID
                     IF (LAC .LT. LBD) THEN
                        LAB = LAC * KEYSQ + LBD
                     ELSE
                        LAB = LBD * KEYSQ + LAC
                     ENDIF
                     WRITE (32+K) JA, npos, LAB, VCOEFF
!xhh - swap not used
!                     IF (K .EQ. 1) WRITE (20) NSWAP
                  ENDIF
               ENDIF
    2       CONTINUE
*
*   All angular coefficients for this pair of CSFs have been
*   generated; update file 30
*
            IF ((JB .EQ. JA) .OR. (.NOT. LINCR))
     :         WRITE (30) JA, JB, npos
    4    CONTINUE
!CFF  if (mod(ja,10) == 0) then
      if (mod(ja,100) == 0) then
!CFF     IF (JA .EQ. 1 .OR. JA .EQ. NCF .OR. MOD (JA,10) .EQ. 0) THEN
         IF (JA .EQ. 1 .OR. JA .EQ. NCF .OR. MOD (JA,100) .EQ. 0) THEN
            PRINT *, 'Row ', JA, '; nnonz = ', npos
     :            ,';  block = ', nb
!         ELSE
!            PRINT *, 'Row ', JA, '; nnonz = ', npos
!     :            ,';  id = ', myid
         ENDIF
      end if
    5 CONTINUE
*
*   Deallocate storage that is no longer required
*
      CALL DALLOC (PNTRIQ)
      CALL DALLOC (PNTJQS)
      CALL DALLOC (PNJCUP)
      CALL ALCBUF (3)
*
*   Write out a report for this run
*
      IF (NDEF .NE. 0 .AND. myid .EQ. 0) THEN
         WRITE (24,*)
         WRITE (24,*) LLISTT, ' T coefficients generated;'
         DO K = 0, KMAX
            WRITE (24,*) LLISTV(K),' V(k=',K,') coefficients generated;'
         ENDDO
      ENDIF
*
*   Set up sparse structure definition arrays in file 30
*
      REWIND (30)

      CALL SETSDA (outsda, npos, LDBPA(4), nb, myid, nprocs, fhead)

      CLOSE (30, STATUS = 'DELETE')
*
*   Sort MCP coefficients into integral-based lists
*
! 1) the T coefficient

      REWIND (31)

      CALL SORT (31, (LLISTT), ntgi, LDBPA(2), nb, fhead)

      IF (NDEF .NE. 0 .AND. myid .EQ. 0)
     &   WRITE (24,*) ntgi, ' I(ab) integrals;'

      CLOSE (31, STATUS = 'DELETE')

! 2) the V coefficient

      DO k = 0, KMAX

         REWIND (k+32)

         CALL SORT (k+32, (LLISTV(k)), ntgi, LDBPA(3), nb, fhead)

         IF (NDEF .NE. 0 .AND. myid .EQ. 0)
     &      WRITE (24,*) ' k = ', k, ': ', ntgi, ' Slater integrals;'

         CLOSE (k+32, STATUS = 'DELETE')

      ENDDO
!xhh - swap not used
!      CLOSE (20, STATUS = 'DELETE')
      CALL DALLOC (PLISTV)

  300 FORMAT (/'From MCP:')
  301 FORMAT (' T_[',1I3,',',1I3,']',
     :   ' (',1I2,1A2,',',1I2,1A2,') = ',1PD19.12)
  302 FORMAT (' V^[(',1I2,')]_[',1I5,',',1I5,']',
     :   ' (',1I2,1A2,',',1I2,1A2,';',
     :        1I2,1A2,',',1I2,1A2,') = ',1PD19.12)

      RETURN
      END
