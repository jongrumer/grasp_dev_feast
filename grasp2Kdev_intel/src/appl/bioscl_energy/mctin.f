************************************************************************
*                                                                      *
      SUBROUTINE MCTIN (IOPAR,JKP,NAME)
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

      include 'parameters.def'
CFF      PARAMETER (NNNW = 120)
CFF      PARAMETER (KEYORB = 121)

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
*
      PARAMETER (CUTOFF = 1.0D-10)
*
      PARAMETER (NFILE = 93)
      PARAMETER (NFILE1 = 237)
*
*   Read the data back as required by OSCL conventions
*
      NFILE2 = NFILE1 +JKP

      M = 0
      K = 0
*
        READ (NFILE2) IBLKI,IBLKF,NW,NKP
        READ (NFILE2) NINT
*
      DO 5 I = 1,NINT
*
           READ (NFILE2) LABL,NCSF
*
         M = M+1
cbieron
         IF (M .GE. NSDIM)
     :      CALL ALCNSA (PNTJJA,PNTJJB,PNTHB1,PNTHB2,
     :                   PNTHC1,PNTHC2,PNTHM1,PNTHM2,
     :                   PNTLAB,PNNPTR,NSDIM,2)
         LAB(M) = LABL
*
*   Read configuration pairs and coefficients for this integral
*
    4    IF (NCSF+K .GT. NTDIM) THEN
            CALL ALCNTA (PISLDR,PISLDR1,PXSLDR,NTDIM,2)
            GOTO 4
         ENDIF
         NPTR(M) = K
           READ (NFILE2) (ISLDR(J+K),ISLDR1(J+K),XSLDR(J+K),J = 1,NCSF)
*          write(*,*) (ISLDR(J+K),XSLDR(J+K),J = 1,NCSF)
         K = K+NCSF
*
    5 CONTINUE
*
*   Close (and hence release) the scratch file
*
c     IF (AVAIL) THEN
c
c        CLOSE (unit=NFILE,status="DELETE")
c      ENDIF
*
      NPTR(M+1) = K
      NINTEG = M
*
      RETURN
*
  301 FORMAT (///1X,I8,' MCT coefficients generated for rank ',I2,
     :        ' and parity ',I2//)
*
      END
