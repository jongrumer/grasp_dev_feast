************************************************************************
*                                                                      *
      SUBROUTINE MATRIX (ncore, j2max)
*       lodcsh2 needs ncore
*                                                                n      *
*   This SUBROUTINE calls routines to  form the  Hamiltonian  matrix   *
*   and to diagonalise it.   The total angular momenta and parity of   *
*   the ASF are found;  the eigenvectors are chosen so that the sign   *
*   of the largest element is positive.                                *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, CONVRT, DALLOC, HMOUT, IQ              *
*                        WGHTD5.                                       *
*               [RCI92]: MANEIG, QED, SETHAM.                          *
*                                                                      *
*                                         Last revision: 28 Dec 1992   *
*   Modified by Xinghong He               Last revision: 31 Dec 1997   *
*   Block version Xinghong He             Last revision: 12 Jun 1998   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)

      integer*8 nelmnt,nelmnt_a


      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
CGG      PARAMETER (NNNWP = 30)
      POINTER (PIENDC,ENDCDUMMY)
      POINTER (PNTRPF,RPFDUMMY)
      POINTER (PNTRQF,RQFDUMMY)
      LOGICAL LFORDR,LTRANS,LVP,LSE,LNMS,LSMS,LDBPG
      CHARACTER*8 CNUM
*
      DIMENSION SLFINT(NNNW) !,UCF(1),SLF_EN(1) COmmented out 2/10/2003
                             ! git and pj
*
      POINTER (PNTSLF,SLF_EN(*))
      POINTER (PNTUCF,UCF(*))

      POINTER (PNETOT,ETOT(*))
*
      POINTER (PNEVAL,EVAL(*))
      POINTER (PNEVEC,EVEC(*))
      POINTER (PNTEMT,EMT(*))
      POINTER (PNIROW,IROW(*))
      POINTER (PNIVEC,IVEC(*))

      POINTER (PNTRIQ,IQA(NNNWP,*))
      POINTER (PNTJQS,JQSA(NNNWP,3,*))
      POINTER (PNJCUP,JCUPA(NNNWP,*))
*
      COMMON/DEBUGG/LDBPG(5)
     :      /DECIDE/LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
     :      /EIGVAL/EAV,PNEVAL
     :      /EIGVEC/PNEVEC
     :      /FOPARM/ICCUT
     :      /HMAT/PNTEMT,PIENDC,PNIROW,NELMNT
     :      /ORB2/NCF,NW,PNTRIQ
     :      /PRNT/NVEC,PNIVEC,NVECMX
     :      /STAT/PNTJQS,PNJCUP
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
     :      /WHERE/IMCDF
     :      /BLIM/IPRERUN,NCSFPRE,COEFFCUT1,COEFFCUT2
     :      /DEF10/AUCM,AUEV,CCMS,FASI,FBSI
     :      /JLABL/JLBL(32),LABJ(32),LABP(2)

*     ...For pre-run
      POINTER (PNEVEC1,EVEC1(*))
      COMMON/EIGVEC1/PNEVEC1

      COMMON/iounit/istdi,istdo,istde

      pointer (pccmin,iccmin(*))
      common/def7/pccmin,notused1,notused2

      POINTER (pncfblk, ncfblk(0:*))
      COMMON/hblock/nblock, pncfblk

      POINTER (pnevblk, nevblk(*))
      POINTER (pncmaxblk, ncmaxblk(*))
      COMMON/hblock2/pnevblk, pncmaxblk

      POINTER (pidxblk, idxblk(*))
      COMMON/blkidx/pidxblk

      POINTER (piccutblk, iccutblk(*))
      COMMON/iccu/piccutblk

      !...nposition+1 is the current position of the .res file
      !...It is set in matrix and used in maneig, spodmv
      COMMON/fposition/nposition

*     ...To deallocate the memory allocated in genintrk.f
      POINTER (PCTEVLRK,VALTEIRK(*))
      POINTER (PCTEILRK, INDTEIRK(*))
      COMMON/CTEILSRK/PCTEILRK,PCTEVLRK

*     ...To deallocate memory allocated in iabint and used in setham
      LOGICAL frstco
      POINTER (PCOEIL,COEILDUMMY)
      POINTER (PCOEVL,COEVLDUMMY)
      COMMON/COEILS/NDCOEA,NCOEI,PCOEIL,PCOEVL,FRSTCO

*     ...To deallocate memory allocated in brint1, brint2,...brint6
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
     :             NDTPA(6),NTPI(6),FIRST(6)  !NTPI,FIRST
      
*     ...To deallocate memory allocated in setham-vpint
      LOGICAL FRSTVP
      POINTER (PINDVP,INDVPDUMMY)
      POINTER (PVALVP,VALVPDUMMY)
      COMMON/VPILST/NDVPA,NVPI,PINDVP,PVALVP,FRSTVP

*     ...To deallocate memory allocated in setham-keint
      LOGICAL FRSTKI
      POINTER (PINDKE,INDKEDUMMY)
      POINTER (PVALKE,VALKEDUMMY)
      COMMON/KEILST/NDKEA,NKEI,PINDKE,PVALKE,FRSTKI

*     ...To deallocate memory allocated in setham-vint
      LOGICAL FRSTVI
      POINTER (PVINIL,VINILDUMMY)
      POINTER (PVINVL,VINVLDUMMY)
      COMMON/VINLST/NDVIN,NVINTI,PVINIL,PVINVL,FRSTVI
!-----------------------------------------------------------------------

      myid = 0
      nprocs = 1
****************************************************************
*     ...Common to all blocks - place here to save CPU time
      CALL auxblk (j2max, atwinv)

****************************************************************
*      Loop over blocks
****************************************************************
      ncminpas = 0

      DO 100 jblock = 1, nblock
         ncf    =   ncfblk(jblock)
         nvec   =   nevblk(jblock)
         nvecmx = ncmaxblk(jblock)
         iccut  = iccutblk(jblock)
         CALL ALLOC (PNTSLF,ncf,8)
         CALL ALLOC (PNTUCF,nw,8)
          do ic=1,ncf
               SLF_EN(IC) = 0.0   
          enddo

         nposition = 7 + nw + nw  ! File position of the previous block
                                  ! in the .res file
         DO i = 1, jblock - 1
            j = (ncfblk(i) - myid - 1 + nprocs) / nprocs
            IF (ncfblk(i) .LT. nprocs) j = ncfblk(i) / (myid+1)
            nposition = nposition + j + 1
         ENDDO

         !.. SETHAM does not handle this extrem case
         IF (nprocs .GT. NCF) 
     &           STOP 'matrix: too many nodes'

*        ...Obtain ivec() from iccmin()
         IF (nvec .GT. 0) THEN
            CALL alloc (pnivec, nvec, 4)
            DO i = 1, nvec
               ivec(i) = iccmin(i+ncminpas)
            ENDDO
            ncminpas = ncminpas + nvec
         ENDIF            

*        ...These 3 were allocated in lodcsh2 and deallocated at the end
*        ... of this routine and in the setham. In this block version, 
*        ... both allocation and deallocation are placed here. See the
*        ... following goto 80 for reason.
         CALL ALLOC (PNTRIQ,NNNWP  *ncf,4)
         CALL ALLOC (PNTJQS,NNNWP*3*ncf,4)
         CALL ALLOC (PNJCUP,NNNWP  *ncf,4)

*	 ...Load CSF list of the current block
			CALL lodcsh2 (21, (ncore), (jblock))
c zou

         IF (LSE) THEN
            PRINT *, 'Entering QED ...'

            WRITE (24,*)
            WRITE (24,*) ' Self Energy Corrections: '
            WRITE (24,*)
               CALL QED_SLFEN (SLFINT)
            DO IC = 1, NCF
               ELEMNT = 0.0D 00
               DO I = 1,NW
c              print *,ic, i,iq(i,ic),slfint(i)
                  ELEMNT = ELEMNT+IQ (I,IC)*SLFINT(I)
               ENDDO
c              print *,ic, elemnt
               SLF_EN(IC) = ELEMNT
            ENDDO

            WRITE (24,*)
            WRITE (24,*) 'Self-energy corrections estimated'
     :              //' --- these will influence the data'
            WRITE (24,*) ' in the RCI92 MIXing coefficients File.'
         ENDIF
c zou

         IF (nvec .LE. 0) THEN
            WRITE (25) jblock, ncf, nvec, 999, 999
*           ...Generate H matrix of the current block and write to disk
*           ...eav, nelmnt are also obtained from genmat
            CALL genmat (atwinv, jblock, myid, nprocs, elsto, irestart,
     :     slf_en)
            	! This call is optional, meaning the matrix of this block
               ! does not need to be computed. But don't comment it out
               ! since other programs may assume thet existence of them.
            CALL genmat2 (irestart, nelmnt_a, elsto)
            GOTO 80	! need to clear memory
         ENDIF
*        ------------------------
         CALL genmat (atwinv, jblock, myid, nprocs, elsto, irestart,
     :     slf_en)
         CALL genmat2 (irestart, nelmnt_a, elsto)
*
*   Allocate and deallocate memory for the mixing coefficients 
*   from the prerun
*
      IF (IPRERUN.EQ.1) CALL ALLOC (PNEVEC1,NCF*NVEC,8)
      IF (IPRERUN.EQ.2) CALL DALLOC (PNEVEC1)

         CALL MANEIG (iiatjpo, iiaspar)
*
*   Write out eigenvalues (ENGOUT), dominant components of the
*   eigenvectors (WGHTD5) to stream 24; write out ASF symmetries,
*   eigenvalues  and eigenvectors to RCI92 mixing coefficients file.
*   EAV and ELSTO are added back to energy here
*
! ELSTO has never been in Hamiltonian matrix, yet it was 
! added to EAV which was later substracted from H. Thus at 
! this point, EAV is correct (it has ELSTO added), EVAL()
! need ELSTO and the correct EAV.
CFF
      IF (NCF > 1) then
         DO i = 1, NVEC
            EVAL(i) = EVAL(i) + ELSTO
         ENDDO
      END IF

      CALL ENGOUT (EAV,EVAL,IiATJPO,iIASPAR,IVEC,NVEC,3)
      CALL WGHTD5 (iiatjpo, iiaspar)

*     ...Write ASF symmetries, eigenvalues, and eigenvectors to RCI92
*     ...MIXing coefficients File; close the file; print a report
      WRITE (25) jblock, ncf, nvec, iiatjpo, iiaspar
      WRITE (25) (ivec(i), i = 1, nvec)
      WRITE (25) EAV,(EVAL(I),I = 1,NVEC)
      WRITE (25) ((EVEC(I+(J-1)*NCF),I = 1,NCF),J = 1,NVEC)
!      print '(5i5)',jblock, ncf, nvec, iiatjpo, iiaspar
!      print '(3I20)',(ivec(i), i = 1, nvec)
!      print '(3f20.15)',EAV,(EVAL(I),I = 1,NVEC)
!      print '(3f20.15)',((EVEC(I+(J-1)*NCF),I = 1,NCF),J = 1,NVEC) 
      PRINT *, 'RCI92 MIXing coefficients File generated.'
*pdbg
*   Save the mixing coefficients from the prerun
*
      IF (IPRERUN .EQ. 1) THEN
         DO J = 1, NVEC
         DO I = 1, NCF
            EVEC1(I+(J-1)*NCF) = EVEC(I+(J-1)*NCF) 
         ENDDO
         ENDDO
      ENDIF
*
*   Estimate diagonal self-energy corrections; the stored
*   eigenvalues and eigenvectors are not modified by these
*   estimates
*
         IF (.not.LSE) THEN
            PRINT *, 'Entering QED ...'
            CALL ALLOC (PNETOT,NVEC,8)

            WRITE (24,*) 
            WRITE (24,*) ' Self Energy Corrections: '
            WRITE (24,*) 
            WRITE (24,301)
            WRITE (24,*) 
  301 FORMAT (' Level  J Parity',7X,'Hartrees',14X,'Kaysers',
     :         16X,'eV' )
  302 FORMAT (1I3,2X,2A4,1P,3D22.14)
            DO J = 1, nvec
               CALL QED (j,SLFINT,UCF)
               ELEMNT = 0.0D 00
               IC = IVEC(J)
               DO I = 1,NW
                  ELEMNT = ELEMNT+UCF(I)*SLFINT(I)
c                 ELEMNT = ELEMNT+IQ (I,IC)*SLFINT(I)
               ENDDO
               ETOT(J) = EVAL(J)+ELEMNT
c
           EAU = ELEMNT  
           ECM = EAU*AUCM
           EEV = EAU*AUEV
           IP = (IIASPAR+3)/2
           WRITE (24,302) j,LABJ(IiATJPO),LABP(IP),EAU,ECM,EEV
c
            ENDDO

            WRITE (24,*)
            WRITE (24,*) 'Self-energy corrections estimated'
     :              //' --- these do not influence the data'
            WRITE (24,*) ' in the RCI92 MIXing coefficients File.'
            CALL ENGOUT (EAV,ETOT,IiATJPO,iIASPAR,IVEC,NVEC,MODE)
c zou       CALL ENGOUT (EAV+elsto,ETOT,IiATJPO,iIASPAR,IVEC,NVEC,MODE)
            CALL dalloc (PNETOT)
         ENDIF

*        ...Locals
         CALL dalloc (pnivec)
*        ...Allocated in maneig
         CALL dalloc (pneval)
         CALL dalloc (pnevec)

   80    CONTINUE

*        ...Locals
         CALL dalloc (PNTRIQ)
         CALL dalloc (PNTJQS)
         CALL dalloc (PNJCUP)
         CALL dalloc (PNTSLF)
         CALL dalloc (PNTUCF)

  100 CONTINUE
*
*   Close the restart files; nothing will be added to them now
*
      CLOSE (imcdf)
*     ...Clean up
      CALL dalloc (pncfblk)
      CALL dalloc (pnevblk)
      CALL dalloc (pncmaxblk)
      CALL dalloc (piccutblk)
      CALL dalloc (pccmin)	! allocated in items as pnccmin

      CALL dalloc (PCTEVLRK)	! allocated in genintrk
      CALL dalloc (PCTEILRK)	! allocated in genintrk
*
*   Deallocate storage for the integral lists from the
*   Dirac-Coulomb operator; the storage was allocated
*   in IABINT and RKINTC
*
      IF (NCOEI .GT. 0) THEN
         CALL DALLOC (PCOEIL)
         CALL DALLOC (PCOEVL)
      ENDIF
*
*   Deallocate storage for the integral lists from the
*   transverse photon interaction operator; this storage
*   was allocated in BRINT1, brint2,...brint6
*
      IF (LTRANS) THEN
         IF (NTPI(1) .GT. 0) THEN
            CALL DALLOC (PINDT1)
            CALL DALLOC (PVALT1)
         ENDIF
         IF (NTPI(2) .GT. 0) THEN
            CALL DALLOC (PINDT2)
            CALL DALLOC (PVALT2)
         ENDIF
         IF (NTPI(3) .GT. 0) THEN
            CALL DALLOC (PINDT3)
            CALL DALLOC (PVALT3)
         ENDIF
         IF (NTPI(4) .GT. 0) THEN
            CALL DALLOC (PINDT4)
            CALL DALLOC (PVALT4)
         ENDIF
         IF (NTPI(5) .GT. 0) THEN
            CALL DALLOC (PINDT5)
            CALL DALLOC (PVALT5)
         ENDIF
         IF (NTPI(6) .GT. 0) THEN
            CALL DALLOC (PINDT6)
            CALL DALLOC (PVALT6)
         ENDIF
      ENDIF
*
*   Deallocate storage for the nuclear motional energy integral
*   lists; this was allocated in KEINT and VINT
*
      IF (LNMS) THEN
         IF (NKEI .GT. 0) THEN
            CALL DALLOC (PINDKE)
            CALL DALLOC (PVALKE)
         ENDIF
      ENDIF
      IF (LSMS) THEN
         IF (NVINTI .GT. 0) THEN
            CALL DALLOC (PVINIL)
            CALL DALLOC (PVINVL)
         ENDIF
      ENDIF
*
*   Deallocate storage for the vacuum-polarisation integral list;
*   this was allocated in VPINT
*
      IF (LVP) THEN
         IF (NVPI .GT. 0) THEN
            CALL DALLOC (PINDVP)
            CALL DALLOC (PVALVP)
         ENDIF
      ENDIF

      CALL dalloc (PNTRPF)	! lodrwf or lodres
      CALL dalloc (PNTRQF)	! lodrwf or lodres

      RETURN

      END
