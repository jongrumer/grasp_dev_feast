************************************************************************
*                                                                      *
      SUBROUTINE BIOTR1(PI,QI,NLI,NINSHLI,
     &                  PF,QF,NLF,NINSHLF,
     &                  NGRID,MXL,SCR,LSCR,NTESTG,
     &                  CISHL,CICI,CFSHL,CFCI)
*                                                                      *
* Generate Matrices for rotating radial functions and                  *
* for counter rotating CI coefficients                                 *
*                                                                      *
* This subroutine is a modification of the corresponding MCHF          *
* routine                                                              *
*                                                                      *
* Per Jonsson,     Lund June 1996                                      *
*                                                                      *
* This subroutine has been modified to support transformations         *
* in cases like 3s2p1d 3s2p                                            *
*                                                                      *
* Per Jonsson,     Lund Feb  1997                                      *
*                                                                      *
* =====                                                                *
* Input                                                                *
* =====                                                                *
*                                                                      *
* PI,QI   : Input radial shells for initial state                      *
*         : (large and small component)                                *
* NLI     : Number of shells per kappa for initial state               *
* NINSHLI : Number of inactive shells  per kappa for initial state     *
*                                                                      *
* PF,QF   : Input radial shells for final   state                      *
*         : (large and small component)                                *
* NLF     : Number of shells per kappa for final   state               *
* NINSHLF : Number of inactive shells  per kappa for final state       *
*                                                                      *
* NGRID   : Number of gridpoints                                       *
* MXL     : Number of kappa quantum numbers                            *
* SCR     :                                                            *
* LSCR    : Total length of scratch space.                             *
* NTESTG  : Global PRINT flag : = 0 => complete silence                *
*                               = 1 => test and PRINT overlap matrices,*
*                                      PRINT header and  t matrix      *
*                            .gt. 1 => Hope you now what you are doing *
*                                                                      *
* ======                                                               *
* Output                                                               *
* ======                                                               *
*                                                                      *
* CISHL : Rotation matrix for initial state shells                     *
* CICI  : Rotation matrix for initial state CI coefficients            *
* CFSHL : Rotation matrix for final state shells                       *
* CFCI  : Rotation matrix for final state CI coefficients              *
* PI,QI : Initial shells in biorthogonal basis                         *
*         biorthogonal basis                                           *
* PF,QF : Final   shells in biorthogonal basis                         *
*                                                                      *
*                                                                      *
* Inactive Shells : The functions of the inactive shells must be       *
*                   be supplied to the program, and NLI,NLF            *
*                   must refer to the total number of                  *
*                   occupied shells ( inactive+active)                 *
*                                                                      *
************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*. Input
*
      DIMENSION NLI(MXL),NLF(MXL)
      DIMENSION NINSHLI(MXL),NINSHLF(MXL)
*
*. Input and output
*
      DIMENSION PI(NGRID,*),QI(NGRID,*),PF(NGRID,*),QF(NGRID,*)
*
*. Output
*
      DIMENSION CISHL(*),CICI(*),CFSHL(*),CFCI(*)
*
*. Scratch
*
      DIMENSION SCR(*)
*
      NTESTL = 0000    
      NTEST = MAX(NTESTL,NTESTG)
      IF(NTEST.GE.1) THEN
        WRITE(6,*)
        WRITE(6,*) '                   *************************'
        WRITE(6,*) '                   *   Entering BIOTR1     *'
        WRITE(6,*) '                   *************************'
        WRITE(6,*)
      END IF
*
*. Scratch should at least be of length 
*
      ILI = 1
      ILF = 1
*
*. Largest number of shells of a given symmetry
*
      NLIMX = IFNMNX(NLI,MXL,1)
      NLFMX = IFNMNX(NLF,MXL,1)
      NLIFMX = MAX(NLIMX,NLFMX)
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' NLIMX,NLFMX NLIFMX ',
     &               NLIMX,NLFMX,NLIFMX
      END IF
*
*. Total numner of shells
*
      NLTI = IELSUM(NLI, MXL)
      NLTF = IELSUM(NLF, MXL)
      IF(NTEST.GE.10)
     &write(6,*) ' NLTI NLTF', NLTI,NLTF
*
* Scratch space for orbital rotations
      KFREE = 1
*
      KLPI = KFREE
      KFREE = KFREE + NLIMX * NGRID
      PRINT*,' In biotrn: KLPI    = ',KLPI
      PRINT*,'            KFREE   = ',KFREE
*
      KLPF = KFREE
      KFREE = KFREE + NLFMX * NGRID
      PRINT*,' In biotrn: KLPF    = ',KLPF
      PRINT*,'            KFREE   = ',KFREE
*
      KLPIF = KFREE
      KFREE = KFREE + NLIFMX * NGRID
      PRINT*,' In biotrn: KLPIF   = ',KLPIF
      PRINT*,'            KFREE   = ',KFREE
*
*. Total overlap matrix
*
      KLSTOT = KFREE
      KFREE = KFREE + NLTI*NLTF
      PRINT*,' In biotrn: KLSTOT  = ',KLSTOT
      PRINT*,'            KFREE   = ',KFREE
*
      KLSIF = KFREE
      KFREE = KFREE + NLIFMX ** 2
      PRINT*,' In biotrn: KLSIF   = ',KLSIF
      PRINT*,'            KFREE   = ',KFREE
*
      KLSIFI = KFREE
      KFREE = KFREE + NLIFMX ** 2
      PRINT*,' In biotrn: KLSIFI  = ',KLSIFI
      PRINT*,'            KFREE   = ',KFREE
*
      KLCI = KFREE
      KFREE = KFREE + NLIFMX ** 2
      PRINT*,' In biotrn: KLCI    = ',KLCI
      PRINT*,'            KFREE   = ',KFREE
*
      KLCF = KFREE
      KFREE = KFREE + NLIFMX ** 2
      PRINT*,' In biotrn: KLCF    = ',KLCF
      PRINT*,'            KFREE   = ',KFREE
*
      KLSCR = KFREE
      KFREE = KFREE + NLIFMX ** 2 + NLIFMX*(NLIFMX+1)
      PRINT*,' In biotrn: KLSCR   = ',KLSCR
      PRINT*,'            KFREE   = ',KFREE
      PRINT*,'         =>  FREE =     ',KFREe
*
*. Check length of scratch
*
      IF(LSCR.LE.KFREE-1) THEN
        WRITE(6,*) ' BIOTR1 in trouble ! '
        WRITE(6,*) 
     &  ' Increase dimension of scratch before call to BIOTR1'
        WRITE(6,*) 
     &  ' Current and required length (LSCR,KFREE-1)',LSCR,KFREE-1
        STOP'Increase LWORK before call to BIOTR1'
      END IF
*
*. Obtain overlap matrix
*
      CALL GETS(SCR(KLSTOT),NLTI,NLTF)

      DO 1000 L = 1, MXL
        IF(NTEST.GE.5) THEN
          WRITE(6,*) '   L = ',L  
          WRITE(6,*) '       Orbital rotation...'
        END IF
*
*. Offsets for given L in shell matrices
*
        IF(L.EQ.1) THEN
          IIOFF = 1
          IFOFF = 1
        ELSE 
          IIOFF = IIOFF + NLI(L-1)** 2
          IFOFF = IFOFF + NLF(L-1)** 2
        END IF
        IF(NTEST.GE.1) WRITE(6,*) 
        IF(NTEST.GE.1) WRITE(6,*) 
     &  ' BIOTRN : Information on transformations of shells with L =',L
        IF(NTEST.GE.1) WRITE(6,*) 
*
* =========================================================
* 1 : Obtain Biorthogonal forms of initial and final shells
* =========================================================
*
*. The overlap matrix can be written
*
*     * * * * * * *
*     *       *   *
*     *   S   * X *
*     *       *   *
*     * * * * * * *
*
* where S is the quadratic subblock.
* The basis functions corresponding to the quadratic
* subblock are made biorthogonal with an UL decomposition.
* The remaining basis
* functions are made biorthogonal by choosing  the
* transformation
*            * * *
*            *   *
*            * Y *
*            *   *
*            * * *
*            * 1 *
*            * * *
*
*
* With Y = -S-1*X
*
        NI = NLI(L)
        NF = NLF(L)
*
*.1.1 : Obtain shells of given L in proper order for
*       biorthogonal treatment
*
        IF(L .EQ. 1 ) THEN
          ILI = 1
          ILF = 1
         ELSE
          ILI = ILI + NLI(L-1)
          ILF = ILF + NLF(L-1)
        END IF

*
* 1.2 obtain biorthogonal of the first min(ni,nf) shells
*
        NIFMN = MIN(NI,NF)
*
*. Overlap matrix SIF = Integral (PI(I)*PF(J))
*
Cww Per change to support cases like 3s2p1d 3s2p. The belonging endif
Cww Is just at the end

       IF (NIFMN.GT.0) THEN

       DO 51 III = 1, NIFMN
       DO 51 JJJ = 1, NIFMN
        SCR(KLSIF+(JJJ-1)*NIFMN+III-1) = 
     &  SCR(KLSTOT-1+(JJJ+ILF-1-1)*NLTI+III+ILI-1)
   51  CONTINUE
       IF(NTEST.GE.15) THEN
         WRITE(6,*) ' Overlap matrix '
         CALL WRTMAT(SCR(KLSIF),NIFMN,NIFMN,NIFMN,NIFMN)
       END IF
*
* Obtain upper triangular CI and CF so CI(T) S CF = 1
* or CF CI(T) = S-1, which corresponds to an UL decomposition
*
*. Invert S
*
         CALL COPVEC(SCR(KLSIF),SCR(KLSIFI),NIFMN**2)
         Call INVMAT(SCR(KLSIFI),SCR(KLCI),NIFMN,NIFMN)

*. UL decompose
         CALL COPVEC(SCR(KLSIFI),SCR(KLSIF),NIFMN**2)
         CALL ULLA(SCR(KLSIF),SCR(KLCF),SCR(KLCI),
     &             NIFMN,SCR(KLSCR))
         CALL TRPMAT(SCR(KLCI),NIFMN,NIFMN,SCR(KLSCR))
         CALL COPVEC(SCR(KLSCR),SCR(KLCI),NIFMN**2)
*
*. The transformation matrix between the first NIFMX
*. shells is now known, biorthogonalize remaining orbitals
*
         IF(NI.NE.NF.AND.NI.NE.0.AND.NF.NE.0) THEN
           IF(NI.GT.NF) THEN
             KLPMX = KLPI
             KLPMN = KLPF
             NMX = NI
             NMN = NF
             KLCMX = KLCI
             KLCMN = KLCF
           ELSE
             KLPMX = KLPF
             KLPMN = KLPI
             NMX = NF
             NMN = NI
             KLCMX = KLCF
             KLCMN = KLCI
           END IF
           NDIFF = NMX - NMN
*
* Y = -S-1 * X
*. overlap X between remaining orbitals and the other set
*
           IF(NI.GT.NF) THEN
*
* I columns F rows
*
             DO 52 III = NMN+1, NMX
             DO 52 JJJ = 1, NF
              SCR(KLSIF+(III-NMN-1)*NF+JJJ-1) = 
     &        SCR(KLSTOT-1+(JJJ+ILF-1-1)*NLTI+III+ILI-1)
   52        CONTINUE
           ELSE IF (NF.GT.NI) THEN
* F columns I rows
             DO 53 JJJ = NMN+1, NMX
             DO 53 III = 1, NI
               SCR(KLSIF+(JJJ-NMN-1)*NI+III-1) = 
     &         SCR(KLSTOT-1+(JJJ+ILF-1-1)*NLTI+III+ILI-1)
   53        CONTINUE
           END IF
*
           IF(NI.GT.NF) THEN
             CALL TRPMAT(SCR(KLSIFI),NMN,NMN,SCR(KLSCR))
             CALL COPVEC(SCR(KLSCR),SCR(KLSIFI),NMN ** 2 )
           END IF
           CALL MATML4(SCR(KLSCR),SCR(KLSIFI),SCR(KLSIF),
     &                 NMN,NDIFF,NMN,NMN,NMN,NDIFF,0)
           CALL SCALVE(SCR(KLSCR),-1.0D0,NMN*NDIFF)
           CALL COPVEC(SCR(KLSCR),SCR(KLSIF),NMN*NDIFF)
*
* Construct complete CMX
*
           CALL SETVEC(SCR(KLSCR),0.0D0,NMX**2)
           DO 300 J = 1, NMX
            IF(J.LE.NIFMN) THEN
             CALL COPVEC(SCR(KLCMX+(J-1)*NIFMN),
     &                   SCR(KLSCR+(J-1)*NMX),NMN)
            ELSE
             CALL COPVEC(SCR(KLSIF+(J-NMN-1)*NMN),
     &                   SCR(KLSCR+(J-1)*NMX),NMN)
             SCR(KLSCR-1+(J-1)*NMX+J) = 1.0D0
            END IF
  300      CONTINUE
*
           CALL COPVEC(SCR(KLSCR),SCR(KLCMX),NMX**2)
         END IF
Cww Pertest
C         ENDIF
*
*. The two upper triangular matrices CI and CF are now known
*. Transfer to permanent arrays
*
         CALL COPVEC(SCR(KLCI),CISHL(IIOFF),NI**2)
         CALL COPVEC(SCR(KLCF),CFSHL(IFOFF),NF**2)
*
*. Rotate the large component of the shells 
*
         CALL COPVEC(PI(1,ILI),SCR(KLPI),NI*NGRID)
         CALL COPVEC(PF(1,ILF),SCR(KLPF),NF*NGRID)
*
         WRITE(*,*) 'Transformation matrices initial'
         CALL WRTMAT(SCR(KLCI),NI,NI,NI,NI)
         CALL MATML4(SCR(KLPIF),SCR(KLPI),SCR(KLCI),
     &               NGRID,NI,NGRID,NI,NI,NI,0)
         CALL COPVEC(SCR(KLPIF),PI(1,ILI),NI*NGRID)
         WRITE(*,*) 'Transformation matrices final'
         CALL WRTMAT(SCR(KLCF),NF,NF,NF,NF)
         CALL MATML4(SCR(KLPIF),SCR(KLPF),SCR(KLCF),
     &               NGRID,NF,NGRID,NF,NF,NF,0)
         CALL COPVEC(SCR(KLPIF),PF(1,ILF),NF*NGRID)
*
*. Rotate the small component of the shells 
*
         CALL COPVEC(QI(1,ILI),SCR(KLPI),NI*NGRID)
         CALL COPVEC(QF(1,ILF),SCR(KLPF),NF*NGRID)
*
         CALL MATML4(SCR(KLPIF),SCR(KLPI),SCR(KLCI),
     &               NGRID,NI,NGRID,NI,NI,NI,0)
         CALL COPVEC(SCR(KLPIF),QI(1,ILI),NI*NGRID)
         CALL MATML4(SCR(KLPIF),SCR(KLPF),SCR(KLCF),
     &               NGRID,NF,NGRID,NF,NF,NF,0)
         CALL COPVEC(SCR(KLPIF),QF(1,ILF),NF*NGRID)
*
         IF(NTEST .GE. 1 ) THEN
           WRITE(6,*) ' Test of overlap of biorthonormal functions'
* F columns I rows
          DO 54 JJJ =1, NF       
          DO 54 III = 1, NI
           SCR(KLSIF+(JJJ-1)*NI+III-1) = 
     &     SCR(KLSTOT-1+(JJJ+ILF-1-1)*NLTI+III+ILI-1)
   54     CONTINUE
          CALL MATML4(SCR(KLSCR),SCR(KLCI),SCR(KLSIF),
     &                NI,NF,NI,NI,NI,NF,1)
          CALL MATML4(SCR(KLSIF),SCR(KLSCR),SCR(KLCF),
     &                NI,NF,NI,NF,NF,NF,0)
           WRITE(6,*)  
     &     ' new overlap matrix ( should be 1 on diag, 0 elsewhere )'
           CALL WRTMAT(SCR(KLSIF),NI,NF,NI,NF)
         END IF

         IF(NTEST.GE.1) THEN
           WRITE(6,*)
           WRITE(6,*) ' Orbital Rotation matrix for I state'
           CALL WRTMAT(CISHL(IIOFF),NI,NI,NI,NI)
           WRITE(6,*) ' Orbital Rotation matrix for F state'
           CALL WRTMAT(CFSHL(IFOFF),NF,NF,NF,NF)
           WRITE(6,*)
         END IF
*
*. Matrix for counterrotation of CI coefficients, initial state
*
         KLTI = KLSIF
         CALL PAMTMT(SCR(KLCI),SCR(KLTI),SCR(KLSCR),NI)
         DO 301 I = 1, NI
           TII = SCR(KLTI-1+(I-1)*NI+I)
           TIII = 1.0D0 / TII
           CALL SCALVE(SCR(KLTI+(I-1)*NI),TIII,I-1)
  301    CONTINUE
         CALL COPVEC(SCR(KLTI),CICI(IIOFF),NI*NI)
*
*. Matrix for counterrotation of CI coefficients, Final state
*
         KLTF = KLSIF
         CALL PAMTMT(SCR(KLCF),SCR(KLTF),SCR(KLSCR),NF)
         DO 302 I = 1, NF
           TII = SCR(KLTF-1+(I-1)*NF+I)
           TIII = 1.0D0 / TII
           CALL SCALVE(SCR(KLTF+(I-1)*NF),TIII,I-1)
  302    CONTINUE
         CALL COPVEC(SCR(KLTF),CFCI(IFOFF),NF*NF)
         IF(NTEST.GE.1) THEN
           WRITE(6,*)
           WRITE(6,*) ' CI-Rotation matrix for I state'
           CALL WRTMAT(CICI(IIOFF),NI,NI,NI,NI)
           WRITE(6,*) ' CI-Rotation matrix for F state'
           CALL WRTMAT(CFCI(IFOFF),NF,NF,NF,NF)
           WRITE(6,*)
         END IF
*
*. End of loop over L
*
Cww Per Change to support case like 3s2p1d 3s2p
         ENDIF

 1000 CONTINUE
*
       RETURN
       END
