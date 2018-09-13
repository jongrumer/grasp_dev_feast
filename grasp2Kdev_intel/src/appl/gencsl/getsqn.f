************************************************************************
*                                                                      *
      SUBROUTINE GETSQN (NCORE,NORB,IOCCS)
*                                                                      *
*   Gets the triads of  subshell quantum numbers  (v, w, 2J)  in the   *
*   arrays  JVLIST, JWLIST, JLIST, each of length  LJL, for all peel   *
*   subshells.                                                         *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 16 Oct 1994   *
*                                                                      *
************************************************************************
*
      PARAMETER (LJLMAX = 20)
*
      include 'parameters.def'
CGG      INTEGER NNNW
CGG      PARAMETER (NNNW = 120)
      DIMENSION IOCCS(NNNW)
*
      COMMON/ORBBOX/JVLIST(NNNW,LJLMAX),JWLIST(NNNW,LJLMAX),
     :              JLIST(NNNW,LJLMAX),LJL(NNNW),IOPTR(NNNW)
     :      /ORBNUM/NP(NNNW),N2J(NNNW),NL(NNNW)
     :      /TERMS/NROWS,ITAB(31),JTAB(32),NTAB(327)
*
      DO 2 IORB = NCORE+1,NORB
*
         IQ = IOCCS(IORB)
         ITJO = N2J(IORB)
*
*   Table in TERMS is compressed
*
         IQT = MIN (IQ,ITJO+1-IQ)
*
*   Read data from appropriate row of TERMS table and load the
*   arrays; set NSS
*
         IF (IQT .GT. 0) THEN
            LOC = (ITJO-1)/2
            LOC = (LOC*(LOC+1))/2+IQT
            NBEG = JTAB(LOC+1)+1
            NEND = JTAB(LOC+2)
            NSS = 0
            DO 1 I = NBEG,NEND,3
               NSS = NSS+1
               JVLIST(IORB,NSS) = NTAB(I  )
               JWLIST(IORB,NSS) = NTAB(I+1)
               JLIST(IORB,NSS)  = NTAB(I+2)-1
    1       CONTINUE
         ELSE
            NSS = 1
            JVLIST(IORB,1) = NTAB(1)
            JWLIST(IORB,1) = NTAB(2)
            JLIST(IORB,1)  = NTAB(3)-1
         ENDIF
*
         LJL(IORB) = NSS
*
    2 CONTINUE
*
      RETURN
      END
