************************************************************************
*                                                                      *
      SUBROUTINE MCTOUT (IOPAR,JKP,NAME)
*                                                                      *
*   This routine loads  coefficients with parity and  rank specified   *
*   by KP(JKP) into the arrays ISLDR and XSLDR.  IOPAR is the parity   *
*   (+/- 1) and is determined from the sign of  KP(JKP).               *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, DALLOC, RALLOC.                        *
*               [OSCL92]: ALCNSA, ALCNTA, TRSORT.                      *
*                                                                      *
*                                         Last revision: 28 Dec 1992   *
*   Updated by Jacek Bieron               Last revision: 10 Mar 1994   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)

      PARAMETER (NNNW = 120)
      PARAMETER (KEYORB = 121)

Cww      INTEGER PNTJJA,PNTJJB,PNTHB1,PNTHB2,PNTHC1,PNTHC2,PNTHM1,PNTHM2,
Cww     :        PNTRIQ
      POINTER (PNTJJA,JJADUMMY)
      POINTER (PNTJJB,JJBDUMMY)
      POINTER (PNTHB1,HB1DUMMY)
      POINTER (PNTHB2,HB2DUMMY)
      POINTER (PNTHC1,HC1DUMMY)
      POINTER (PNTHC2,HC2DUMMY)
      POINTER (PNTHM1,HM1DUMMY)
      POINTER (PNTHM2,HM2DUMMY)
      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL LDBPA,LFORDR,LTRANS,LVP,LSE,LNMS,LSMS,AVAIL
      CHARACTER*2 NH
      CHARACTER*24 NAME(2)
*
      DIMENSION TSHELL(NNNW)
*
      POINTER (PLABEL,LABEL(1)),(PCOEFF,COEFF(1))
*
      POINTER (PISLDR,ISLDR(1)),(PXSLDR,XSLDR(1))
      POINTER (PISLDR1,ISLDR1(1))
      POINTER (PNTLAB,LAB(1)),(PNNPTR,NPTR(1))
      POINTER (PNTRKP,KP(1))
*
      COMMON/DEBUGA/LDBPA(5)
     :      /DECIDE/LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
     :      /FOPARM/ICCUT
     :      /OFFD/NOFFD1,NOFFD2
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /OSC1/PNTJJA,PNTJJB,PNTHB1,PNTHB2,
     :            PNTHC1,PNTHC2,PNTHM1,PNTHM2,NSDIM
     :      /OSC2/LK,KK
     :      /OSC3/PXSLDR,PISLDR,PISLDR1,NTDIM
     :      /OSC5/NINT,PNTLAB,PNNPTR,NINTEG
     :      /OSC6/NKP,PNTRKP

      include 'mpif.h'
      integer  myid, nprocs, ierr
      COMMON /mpi/ myid, nprocs, ierr                                           
*
      PARAMETER (CUTOFF = 1.0D-10)
*
      PARAMETER (NFILE = 93)
      PARAMETER (NFILE1 = 237)
*
* NCFI(I): the end position of the Ith block for the initial states in the globle CSF list
* NCFF(I): the end position of the Ith block for the final states in the globle CSF list
CFF      common /IBLK/NBLOCKI,NCFI(10)
CFF      common /FBLK/NBLOCKF,NCFF(10)
      common /IBLK/NBLOCKI,NCFI(50)
      common /FBLK/NBLOCKF,NCFF(50)
*
*   Check if angular data is available on file
*
      NFILE2 = NFILE1 +JKP
      CALL ANGDATA(NAME,AVAIL,JKP,NFILE2)
*
*   If angular data is not available open the scratch file to store the MCT
*   coefficients; position file
*   to beginning
*
*
      LK = ABS (KP(JKP))
      IOPAR = ISIGN (1,KP(JKP))
      IF (AVAIL) RETURN
*
* Start of the block loops
      do 1001 iblki = 1, nblocki
      do 1000 iblkf = 1, nblockf
*
*       OPEN (NFILE,STATUS = 'new', FORM = 'UNFORMATTED')
        OPEN (NFILE,STATUS = 'unknown', FORM = 'UNFORMATTED')
        NMCT = 0
*
*   
*   If angular data is not available
*   Generate MCT coefficients for given rank/parity combination and
*   store them by CSF on NFILE
*
*
*   Allocate storage to buffer arrays
*
        NLABEL = 32
        CALL ALLOC (PLABEL,NLABEL,4)
        CALL ALLOC (PCOEFF,NLABEL,8)

        IF(iblki.eq.1) then
          NCFI0 = 1
        else
          NCFI0 = NCFI(iblki-1)+1
        endif
c
        IF(iblkf.eq.1) then
          NCFF0 = NCFI(NBLOCKI)+1
        else
          NCFF0 = NCFI(NBLOCKI)+NCFF(iblkf-1)+1
        endif
        DO 3 IC = myid+NCFI0, NCFI(IBLKI), nprocs
        DO 3 IR = NCFF0,NCFI(NBLOCKI)+NCFF(iblkf)
*
*          IR = IC
*
           NCR = 0

*                                                                        
*   In many case one is interested only in M1 and E2 transitions between
*   levels with different J values. If this is the case then the do check
*   on the J quantum numbers of the CSFs before calling TNSRJJ.         
*                                                                       
           IF (KP(JKP).EQ.1.AND.NOFFD1.EQ.1) THEN                         
             IF (ITJPO (IC).EQ.ITJPO (IR)) GOTO 13                        
           ENDIF                                                          
           IF (KP(JKP).EQ.2.AND.NOFFD2.EQ.1) THEN                         
             IF (ITJPO (IC).EQ.ITJPO (IR)) GOTO 13                        
           ENDIF
c          if(ispar(ic)*ispar(ir)*iopar.ne.1.
c    &        or.itrig(itjpo(ic),itjpo(ir),2*lk+1).ne.1) go to 13
c          if(ichkq1(IC,IR).eq.0) go to 13

           CALL TNSRJJ (LK,IOPAR,IC,IR,IA,IB,TSHELL)
           IF (IA .NE. 0) THEN
              IF (IA .EQ. IB) THEN
                 DO 2 IA = 1,NW
                    IF (ABS (TSHELL(IA)) .GT. CUTOFF) THEN
                       NCR = NCR+1
                       IF (NCR .GT. NLABEL) THEN
                          NEWSIZ = 2*NLABEL
                          CALL RALLOC (PLABEL,NLABEL,NEWSIZ,4)
                          CALL RALLOC (PCOEFF,NLABEL,NEWSIZ,8)
                          NLABEL = NEWSIZ
                       ENDIF
                       LABEL(NCR) = IA*KEYORB+IA
                       COEFF(NCR) = TSHELL(IA)
                    ENDIF
    2            CONTINUE
              ELSE
                 IF (ABS (TSHELL(1)) .GT. CUTOFF) THEN
                    NCR = NCR+1
                    IF (NCR .GT. NLABEL) THEN
                       NEWSIZ = 2*NLABEL
                       CALL RALLOC (PLABEL,NLABEL,NEWSIZ,4)
                       CALL RALLOC (PCOEFF,NLABEL,NEWSIZ,8)
                       NLABEL = NEWSIZ
                    ENDIF
                    LABEL(NCR) = IA*KEYORB+IB
                    COEFF(NCR) = TSHELL(1)
                 ENDIF
              ENDIF
           ENDIF
           IF (NCR .GT. 0) THEN
              WRITE (NFILE) IC-NCFI0+1,IR-NCFF0+1,NCR
c             WRITE (NFILE) IC,IR,NCR
              WRITE (NFILE) (LABEL(I),COEFF(I),I = 1,NCR)
              NMCT = NMCT+NCR
           ENDIF
*
*                                                                       
   13      CONTINUE                                                       
                                                                        
*
    3   CONTINUE
*
*   Deallocate storage for buffer arrays
*
        CALL DALLOC (PLABEL)
        CALL DALLOC (PCOEFF)
*
C       WRITE (*,301) NMCT,LK,IOPAR
*
*   Sort the MCT coefficients by integral labels
*
        CALL TRSORT (NAME,NFILE,NFILE2,LDBPA(2),JKP,IBLKI,IBLKF)
        close (NFILE,status='delete')
* end of the loops for blocks
 1000 continue
 1001 continue
*
*   Read the data back as required by OSCL conventions
*
        REWIND (NFILE2)
      return
*
  301 FORMAT (///1X,I8,' MCT coefficients generated for rank ',I2,
     :        ' and parity ',I2//)
*
      END
