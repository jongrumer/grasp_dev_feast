************************************************************************
*                                                                      *
      SUBROUTINE GENNRL
*                                                                      *
*   Generates the CSF list using nonrelativistic subshell labels.      *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, DALLOC, RALLOC, RALC2D.                *
*               [GENCSL]: DSTBUT, MAKCSR, PARSJL, PRNTCN, PRNTPJ,      *
*                         PRSNCL, PRSNSL, WRTCSF.                      *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 16 Oct 1994   *
*                                                                      *
************************************************************************
*
!$Id: gennrl.f,v 1.1 2003/09/30 05:57:04 georgio Exp $
!$Log: gennrl.f,v $
!Revision 1.1  2003/09/30 05:57:04  georgio
!
!added
!
!Revision 1.2  1997/06/02 22:02:55  xhh
!*** empty log message ***
!
! Short output lines joined
! XHH 1997.01.21

      LOGICAL DOCHEK,GETYN,PRNTSC,YES
      CHARACTER*2 SYM
      CHARACTER*1 NRSYM
*
*   The maximum number of states of a subshell is LJLMAX; this
*   value is set to the largest value of any element of the
*   array ITAB in BLOCK DATA TERM
*
      PARAMETER (LJLMAX = 20)
*
      include 'parameters.def'
CGG      INTEGER NNNW
CGG      PARAMETER (NNNW = 120)
CGG      INTEGER NNNTJV
CGG      PARAMETER (NNNTJV = 10)
CGG      INTEGER NNNWM1
CGG      PARAMETER (NNNWM1 = 119)
CGG      INTEGER NNNWM2
CGG      PARAMETER (NNNWM2 = 118)
*
      POINTER (PIOREF,IOREF(NNNW,*))
      POINTER (PNIPAR,IPAR(*))
      POINTER (PNTITJ,ITJ(NNNTJV,*))
      POINTER (PNTNTJ,NTJ(*))
*
      POINTER (PJCLST,JCLIST(NNNWM1,*))
*
      DIMENSION IOCNR(NNNW),MNOCNR(NNNW),MXOCNR(NNNW)
      DIMENSION IOCCS(NNNW),MXOCC(NNNW)
      DIMENSION ITJREF(NNNTJV)
*
      COMMON/COUBOX/PJCLST,MNJVC,LJCL(NNNWM1),ICPTR(NNNWM2)
     :      /NRORBN/NPNR(NNNW),N2LNR(NNNW)
     :      /NRORBS/NRSYM(NNNW)
     :      /ORBBOX/JVLIST(NNNW,LJLMAX),JWLIST(NNNW,LJLMAX),
     :              JLIST(NNNW,LJLMAX),LJL(NNNW),IOPTR(NNNW)
     :      /ORBNUM/NP(NNNW),N2J(NNNW),NL(NNNW)
     :      /ORBSYM/SYM(NNNW)
      COMMON/iounit/istdi,istdo,istde
*
    1 WRITE(istde,*) 'Enter the list of nonrelativistic core subshells:'
      NORBNR = 0
      CALL PRSNSL (NORBNR)
      NCOREN = NORBNR
*
      WRITE(istde,*) 'Enter the list of nonrelativistic peel subshells:'
      CALL PRSNSL (NORBNR)
*
*   Ensure that there are no repetitions in the subshell list
*
      DO 3 I = 1,NORBNR-1
         NPREF = NPNR(I)
         N2LREF = N2LNR(I)
         DO 2 J = I+1,NORBNR
            IF ((NPNR(J) .EQ. NPREF) .AND.
     :          (N2LNR(J) .EQ. N2LREF)) THEN
                WRITE(istde,*) 'GENNRL: Detected a repeated subshell;'
     &, ' redo ...'
               GOTO 1
            ENDIF
    2    CONTINUE
    3 CONTINUE
*
*   Determine the maximum nonrelativistic occupation of each subshell;
*   override incorrect user input
*
      IF (NORBNR .GT. NCOREN) THEN
         WRITE(istde,*) 'Redefine default maximum nonrelativistic'
     &, ' peel subshell occupations?'
         YES = GETYN ()
      ELSE
         YES = .FALSE.
      ENDIF
*
      IF (YES) THEN
         WRITE(istde,*) 'Enter the list of maximum occupations'
     &, ' to be redefined:'
         CALL PRSNCN (NPNR,NRSYM,NCOREN,NORBNR,IOCNR)
      ELSE
         DO 4 I = 1,NORBNR
            IOCNR(I) = 0
    4    CONTINUE
      ENDIF
*
*   Store the maximum nonrelativistic occupation;
*   generate the relativistic subshell data
*
      NORB = 0
      DO 5 I = 1,NORBNR
         IF (IOCNR(I) .GT. 0) THEN
            IF (N2LNR(I) .LE. 6) THEN
               MXOCNR(I) = MIN (IOCNR(I),2*(N2LNR(I)+1))
            ELSEIF (N2LNR(I) .EQ. 8) THEN
               MXOCNR(I) = MIN (IOCNR(I),10)
            ELSE
               MXOCNR(I) = MIN (IOCNR(I),4)
            ENDIF
         ELSE
            IF     (N2LNR(I) .LE. 6) THEN
               MXOCNR(I) = 2*(N2LNR(I)+1)
            ELSEIF (N2LNR(I) .EQ. 8) THEN
               MXOCNR(I) = 10
            ELSE
               MXOCNR(I) = 4
            ENDIF
         ENDIF
         IF (NRSYM(I) .NE. 's') THEN
            NORB = NORB+1
            NP(NORB) = NPNR(I)
            N2J(NORB) = N2LNR(I)-1
            NL(NORB) = N2LNR(I)/2
            SYM(NORB) = NRSYM(I)//'-'
            IF (N2J(NORB) .LE. 7) THEN
               MXOCC(NORB) = MIN (N2J(NORB)+1,MXOCNR(I))
            ELSE
               MXOCC(NORB) = MIN (2,MXOCNR(I))
            ENDIF
            NORB = NORB+1
            NP(NORB) = NPNR(I)
            N2J(NORB) = N2LNR(I)+1
            NL(NORB) = N2LNR(I)/2
            SYM(NORB) = NRSYM(I)//' '
            IF (N2J(NORB) .LE. 7) THEN
               MXOCC(NORB) = MIN (N2J(NORB)+1,MXOCNR(I))
            ELSE
               MXOCC(NORB) = MIN (2,MXOCNR(I))
            ENDIF
         ELSE
            NORB = NORB+1
            NP(NORB) = NPNR(I)
            N2J(NORB) = 1
            NL(NORB) = 0
            SYM(NORB) = 's '
            MXOCC(NORB) = MIN (2,MXOCNR(I))
         ENDIF
         IF (I .EQ. NCOREN) NCORE = NORB
    5 CONTINUE
      NORBM1 = NORB-1
      NORBM2 = NORB-2
*
      IF (YES) THEN
         WRITE(istde,*) 'Equivalent maximum relativistic peel'
     &, ' subshell occupations:'
         CALL PRNTCN (NP,SYM,MXOCC,NCORE,NORB)
      ENDIF
*
*   Determine the minimum nonrelativistic occupation of each subshell
*
      IF (NORBNR .GT. NCOREN) THEN
         WRITE(istde,*) 'Redefine default minimum nonrelativistic'
     &, ' peel subshell occupations?'
         YES = GETYN ()
      ELSE
         YES = .FALSE.
      ENDIF
*
      IF (YES) THEN
        WRITE(istde,*) 'Enter list of minimum occupations'
     &, ' to be redefined:'
        CALL PRSNCN (NPNR,NRSYM,NCOREN,NORBNR,MNOCNR)
      ELSE
        DO 6 I = 1,NORBNR
           MNOCNR(I) = 0
    6   CONTINUE
      ENDIF
*
*   Read the reference configurations; nonrelativistic
*   configurations are read, and all equivalent relativistic
*   configurations are generated; determine the number of peel
*   electrons and the atomic parities and J values
*
      NREFD = 2
      CALL ALLOC (PIOREF,NNNW*NREFD,4)
      CALL ALLOC (PNIPAR,       NREFD,4)
      CALL ALLOC (PNTITJ,NNNTJV*NREFD,4)
      CALL ALLOC (PNTNTJ,       NREFD,4)
      NREF = 0
      IF (NORBNR .GT. NCOREN) THEN
    7    WRITE(istde,*) 'Enter a nonrelativistic reference peel'
     &, ' subshell configuration (null if done):'
         CALL PRSNCN (NPNR,NRSYM,NCOREN,NORBNR,IOCNR)
         DO 8 I = NCOREN+1,NORBNR
            IF (IOCNR(I) .GT. 0) GOTO 9
    8    CONTINUE
         GOTO 17
    9    LSUM = 0
         ICOUNT = 0
         DO 10 I = NCOREN+1,NORBNR
            IOCCI = IOCNR(I)
            IF ((IOCCI .GT. MXOCNR(I)) .OR.
     :          (IOCCI .LT. MNOCNR(I))) THEN
               WRITE(istde,*) 'GENNRL: Subshell occupation out of'
     &, ' bounds; redo ...'
               GOTO 7
            ENDIF
            LSUM = LSUM+IOCCI*N2LNR(I)/2
            ICOUNT = ICOUNT+IOCCI
            IPTY = (-1)**LSUM
   10    CONTINUE
         IF (ICOUNT .GT. 0) THEN
            IF (NREF .EQ. 0) THEN
               NPEEL = ICOUNT
            ELSE
               IF (ICOUNT .NE. NPEEL) THEN
                  WRITE(istde,*) 'GENNRL: Inconsistent number of'
     &                   ,' peel electrons. redo ...'
                  GOTO 7
               ENDIF
            ENDIF
            IF (IPTY .EQ. -1) THEN
               WRITE(istde,*) 'Parity of configuration is odd'
            ELSE
               WRITE(istde,*) 'Parity of configuration is even'
            ENDIF
            WRITE(istde,*) 'Enter the J values for this configuration'
     &, ' (null for all possible J values):'
            CALL PARSJL (NPEEL,ITJREF,NTJREF)
            LOC = NCORE
            LVL = NORB
            DO 12 I = NCOREN+1,NORBNR
               IOCNRI = IOCNR(I)
               IF     (N2LNR(I) .EQ. 0) THEN
                  LOC = LOC+1
                  JWLIST(LOC,1) = IOCNRI
                  LJL(LOC) = 1
                  IOPTR(LOC) = 1
               ELSE
                  LOC = LOC+2
                  IOCM = MXOCC(LOC-1)
                  IOCP = MXOCC(LOC  )
                  IF (IOCNRI .LE. IOCM) THEN
                     IRANGE = IOCNRI+1
                     IOCLO = IOCNRI
                  ELSE
                     IRANGE = MIN (IOCP+IOCM-IOCNRI,IOCM)+1
                     IOCLO = IOCM
                  ENDIF
                  DO 11 II = 1,IRANGE
                     JWLIST(LOC-1,II) = IOCLO
                     JWLIST(LOC  ,II) = IOCNRI-IOCLO
                     IOCLO = IOCLO-1
   11             CONTINUE
                  LJL(LOC-1) = IRANGE
                  LJL(LOC  ) = IRANGE
                  IOPTR(LOC-1) = 1
                  IOPTR(LOC  ) = 1
               ENDIF
   12       CONTINUE
   13       NREF = NREF+1
            IF (NREF .GT. NREFD) THEN
               NEWSIZ = NREF+NREF/2
               CALL RALC2D (PIOREF,NNNW,NREFD,NNNW,NEWSIZ,4)
               CALL RALLOC (PNIPAR,       NREFD,       NEWSIZ,4)
               CALL RALC2D (PNTITJ,NNNTJV,NREFD,NNNTJV,NEWSIZ,4)
               CALL RALLOC (PNTNTJ,       NREFD,       NEWSIZ,4)
               NREFD = NEWSIZ
            ENDIF
            DO 14 I = NCORE+1,NORB
               IOREF(I,NREF) = JWLIST(I,IOPTR(I))
   14       CONTINUE
            IPAR(NREF) = IPTY
            NTJ(NREF) = NTJREF
            DO 15 I = 1,NTJREF
               ITJ(I,NREF) = ITJREF(I)
   15       CONTINUE
   16       IF (IOPTR(LVL) .LT. LJL(LVL)) THEN
               IF (N2J(LVL) .EQ. 1) THEN
                  IOPTR(LVL) = IOPTR(LVL)+1
               ELSE
                  IOPTR(LVL  ) = IOPTR(LVL  )+1
                  IOPTR(LVL-1) = IOPTR(LVL-1)+1
               ENDIF
               LVL = NORB
               GOTO 13
            ELSE
               IF (LVL .GT. NCORE) THEN
                  IF (N2J(LVL) .EQ. 1) THEN
                     LVL = LVL-1
                     IOPTR(LVL+1) = 1
                  ELSE
                     LVL = LVL-2
                     IOPTR(LVL+1) = 1
                     IOPTR(LVL+2) = 1
                  ENDIF
                  GOTO 16
               ENDIF
            ENDIF
            GOTO 7
         ENDIF
*
   17    CONTINUE 

Cww   17    WRITE(istde,*) 'Copy list of equivalent relativistic'
Cww         WRITE(istde,*) ' reference configurations to screen?'
Cww         YES = GETYN ()
Cww         IF (YES) THEN
Cww            DO 18 I = 1,NREF
Cww               CALL PRNTCN (NP,SYM,IOREF(1,I),NCORE,NORB)
Cww               CALL PRNTPJ (IPAR(I),NTJ(I),ITJ(1,I))
Cww   18       CONTINUE
Cww         ENDIF
*
*   Determine the maximum number of substitutions allowed
*
   19    WRITE(istde,*) 'Enter the maximum number of substitutions'
     &, ' permitted (0, ..., 9):'
         READ (*,'(1I1)',IOSTAT = IOS) NSUB
         IF ((IOS .NE. 0) .OR. (NSUB .LT. 0)) THEN
            WRITE(istde,*) 'GENRL: Unacceptable response;'
     &, ' redo ...'
            GOTO 19
         ENDIF
*
*   Is a copy of the CSF list to appear on screen?
*
Cww         WRITE(istde,*) 'Copy CSF list to screen?'
Cww         YES = GETYN ()
Cww         IF (YES) THEN
Cww            PRNTSC = .TRUE.
Cww         ELSE
Cww            PRNTSC = .FALSE.
Cww         ENDIF
*
      ENDIF
*
*   Write the header of grasp92.csf
*
      REWIND (21)
      WRITE (21,'(A)') 'Core subshells:'
      WRITE (21,301) (NP(I),SYM(I),I = 1,NCORE)
      WRITE (21,'(A)') 'Peel subshells:'
      IF ((NORB-NCORE).LE.80) THEN
         WRITE (21,301) (NP(I),SYM(I),I = NCORE+1,NORB)
      ELSE
         WRITE (21,301) (NP(I),SYM(I),I = NCORE+1,NCORE+80)
         WRITE (21,301) (NP(I),SYM(I),I = NCORE+81,NORB)
      ENDIF
      WRITE (21,'(A)') 'CSF(s):'
*
*   Determine the location of the last nonempty subshell
*   in the reference distributions
*
      DO 21 I = NCORE+1,NORB
         DO 20 II = 1,NREF
            IF (IOREF(I,II) .NE. 0) ICHECK = I
   20    CONTINUE
   21 CONTINUE
      ICHECK = ICHECK+1
      DOCHEK = ICHECK .LT. NORB
*
*   Initial allocation for array JCLIST
*
      MNJVC = 2
      CALL ALLOC (PJCLST,NNNWM1*MNJVC,4)
*
*   Generate the CSF list
*
      ICOUNT = 0
*
*   Initialise the electron distribution
*
      DO 22 I = 1,NCORE
         IOCCS(I) = MXOCC(I)
   22 CONTINUE
*
      NLEFT = NPEEL
      DO 23 I = NCORE+1,NORB
         IF (NLEFT .GT. 0) THEN
            IOCCS(I) = MIN (MXOCC(I),NLEFT)
            NLEFT = NLEFT-IOCCS(I)
            NLAST = I
         ELSE
            IOCCS(I) = 0
         ENDIF
   23 CONTINUE
*
      IF (NPEEL .GT. 0) THEN
         GOTO 27
      ELSE
         ICOUNT = 1
         GOTO 43
      ENDIF
*
   24 IF (IOCCS(NORB) .EQ. NPEEL) GOTO 43
*
*   Distribute the electrons in the peel subshells
*
      CALL DSTBUT (IOCCS,NCORE,NORB,NLAST)
*
      IF (DOCHEK) THEN
*
*   Expedite the search by omitting all distributions with
*   more than NSUB electrons outside the last subshell
*   occupied in the reference distributions
*
         IF (IOCCS(ICHECK) .GT. NSUB) THEN
            IOCCS(NORB) = IOCCS(ICHECK)
            IOCCS(ICHECK) = 0
            NLAST = NORB
            GOTO 24
         ENDIF
*
      ENDIF
*
*   Reject this distribution if any peel subshell occupation
*   is excessive or insufficient
*
   27 LOC = NCORE
      DO 28 I = NCOREN+1,NORBNR
         IF (N2LNR(I) .GT. 0) THEN
            LOC = LOC+2
            IOCM = IOCCS(LOC-1)
            IOCP = IOCCS(LOC  )
            IOCT = IOCM+IOCP
            IF ((IOCM .GT. MXOCC(LOC-1)) .OR.
     :          (IOCP .GT. MXOCC(LOC  )) .OR.
     :          (IOCT .GT. MXOCNR(I)) .OR.
     :          (IOCT .LT. MNOCNR(I))) GOTO 24
         ELSE
            LOC = LOC+1
            IOCT = IOCCS(LOC)
            IF ((IOCT .GT. MXOCNR(I)) .OR.
     :          (IOCT .LT. MNOCNR(I))) GOTO 24
         ENDIF
   28 CONTINUE
*
*   Compute the parity of this distribution; the parity of the core
*   is always even
*
      LSUM = 0
      DO 29 I = NCORE+1,NORB
         LSUM = LSUM+IOCCS(I)*NL(I)
   29 CONTINUE
      IPTY = (-1)**LSUM
*
*   Reject this distribution if the number of excitations or
*   parity is incorrect with respect to all reference
*   distributions
*
      DO 31 II = 1,NREF
         NEX = 0
         DO 30 I = NCORE+1,NORB
            NEX = NEX + MAX (0,IOREF(I,II)-IOCCS(I))
   30    CONTINUE
         IF (NEX .LE. NSUB) THEN
            IF (IPAR(II) .EQ. IPTY) GOTO 32
         ENDIF
   31 CONTINUE
      GOTO 24
*
*   Determine all peel subshell quantum numbers
*
   32 CALL GETSQN (NCORE,NORB,IOCCS)
*
*   Reset the pointers for the lists of subshell quantum numbers
*
      DO 33 I = NCORE+1,NORB
         IOPTR(I) = 1
   33 CONTINUE
*
*   Determine the first set of all intermediate angular momenta
*
      DO 34 I = 1,NORBM2-NCORE
         ICPTR(I) = 1
         CALL COUPLE (NCORE,I)
   34 CONTINUE
*
*   Initialize level indicator
*
      LVL = NORBM1-NCORE
*
*   Generate all possible final angular momenta
*
   35 DO 36 I = LVL,NORBM1-NCORE
         CALL COUPLE (NCORE,I)
   36 CONTINUE
*
*   Reset level indicator
*
      LVL = NORBM1-NCORE
*
      DO 41 I = 1,LJCL(NORBM1-NCORE)
*
*   Store information only for cases that have the correct J/parity
*   and the correct number of excitations
*
         JATOM = JCLIST(NORBM1-NCORE,I)
         DO 39 II = 1,NREF
            NEX = 0
            DO 37 IV = NCORE+1,NORB
               NEX = NEX + MAX (0,IOREF(IV,II)-IOCCS(IV))
   37       CONTINUE
            IF (NEX .LE. NSUB) THEN
               IF (IPAR(II) .EQ. IPTY) THEN
                  IF (NTJ(II) .EQ. -1) GOTO 40
                  DO 38 III = 1,NTJ(II)
                     IF (ITJ(III,II) .EQ. JATOM) GOTO 40
   38             CONTINUE
               ENDIF
            ENDIF
   39    CONTINUE
*
*   Match not found; skip to end of loop
*
         GOTO 41
*
*   Match found
*
   40    ICOUNT = ICOUNT+1
         CALL WRTCSF (NCORE,NORB,IOCCS,JATOM,IPTY,PRNTSC)
*
   41 CONTINUE
*
*   Reset the pointers and determine the level indicator
*
   42 LVLP1 = LVL+1
      LVLM1 = LVL-1
*
*   Reset the pointer for the current outermost subshell
*
      IOPTR(NCORE+LVLP1) = IOPTR(NCORE+LVLP1)+1
*
      IF (IOPTR(NCORE+LVLP1) .GT. LJL(NCORE+LVLP1)) THEN
*
*   All J values for the outermost subshell have been exhausted;
*   reset the subshell pointer
*
         IOPTR(NCORE+LVLP1) = 1
*
         IF (LVLM1 .GE. 1) THEN
*
*   At this level an intermediate angular momentum is coupled to
*   a subshell angular momentum; increment the pointer for the
*   subshell angular momentum
*
            ICPTR(LVLM1) = ICPTR(LVLM1) + 1
*
            IF (ICPTR(LVLM1) .GT. LJCL(LVLM1)) THEN
*
*   All intermediate J values have been exhausted; reset the pointer
*   for the intermediate angular momentum; go back a level
*
               ICPTR(LVLM1) = 1
               LVL = LVLM1
               GOTO 42
*
            ENDIF
*
         ELSE
*
*   At this level the first two subshells are coupled
*
            IOPTR(NCORE+1) = IOPTR(NCORE+1)+1
*
*   All couplings for this distribution have been exhausted;
*   start with a new distribution
*
           IF (IOPTR(NCORE+1) .GT. LJL(NCORE+1)) GOTO 24
*
         ENDIF
      ENDIF
*
*   Obtain a new set of couplings with the same distribution
*
      GOTO 35
*
*   Print summary
*
   43 WRITE(istde,*) NORB,' relativistic subshells. '
     &, ICOUNT,' relativistic CSFs generated.'
*
*   All done
*
      CALL DALLOC (PIOREF)
      CALL DALLOC (PNIPAR)
      CALL DALLOC (PNTITJ)
      CALL DALLOC (PNTNTJ)
*
      RETURN
*
  300 FORMAT (32(1X,1I2,1A1,'(',1I2,')'))
  301 FORMAT (120(1X,1I2,1A2))
*
      END
