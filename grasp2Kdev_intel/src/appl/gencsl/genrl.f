************************************************************************
*                                                                      *
      SUBROUTINE GENRL
*                                                                      *
*   Generates the CSF list using relativistic subshell labels.         *
*                                                                      *
*   Subprogram(s) called: ALLOC, COUPLE, DALLOC, DSTBUT, PARSJL,       *
*                         PRNTCN, PRNTPJ, PRSRCN, PRSRSL, RALC2D,      *
*                         RALLOC, WRTCSF.                              *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 25 Sep 1993   *
*                                                                      *
************************************************************************
*

!$Id: genrl.f,v 1.1 2003/09/30 05:57:04 georgio Exp $
!$Log: genrl.f,v $
!Revision 1.1  2003/09/30 05:57:04  georgio
!
!added
!
!Revision 1.2  1997/06/02 22:02:55  xhh
!*** empty log message ***
!
      CHARACTER*(*)    RCSID
      PARAMETER        ( RCSID
     & ='$Id: genrl.f,v 1.1 2003/09/30 05:57:04 georgio Exp $'
     & )

! RCSID etc. added
! Short output lines joined
! XHH 1997.01.21

      LOGICAL GETYN,PRNTSC,YES
      CHARACTER*2 SYM
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
      DIMENSION IOCCS(NNNW),MNOCC(NNNW),MXOCC(NNNW)
*
      COMMON/COUBOX/PJCLST,MNJVC,LJCL(NNNWM1),ICPTR(NNNWM2)
     :      /ORBBOX/JVLIST(NNNW,LJLMAX),JWLIST(NNNW,LJLMAX),
     :              JLIST(NNNW,LJLMAX),LJL(NNNW),IOPTR(NNNW)
     :      /ORBNUM/NP(NNNW),N2J(NNNW),NL(NNNW)
     :      /ORBSYM/SYM(NNNW)
       COMMON/iounit/istdi,istdo,istde
*
*   Determine the core and peel subshells
*
    1 WRITE(istde,*) 'Enter the list of relativistic core'
     &, ' subshells'
      NORB = 0
      CALL PRSRSL (NORB)
*
      NCORE = NORB
*
      WRITE(istde,*) 'Enter the list of relativistic peel'
     &, ' subshells:'
      CALL PRSRSL (NORB)
*
      NORBM1 = NORB-1
      NORBM2 = NORB-2
*
*   Ensure that there are no repetitions in the subshell list
*
      DO 3 I = 1,NORB-1
         NPREF = NP(I)
         N2JREF = N2J(I)
         NLREF = NL(I)
         DO 2 J = I+1,NORB
            IF ((NP(J) .EQ. NPREF) .AND.
     :          (N2J(J) .EQ. N2JREF) .AND.
     :          (NL(J) .EQ. NLREF)) THEN
               WRITE(istde,*) 'GENRL: Detected a repeated subshell;'
     &, ' redo ...'
               GOTO 1
            ENDIF
    2    CONTINUE
    3 CONTINUE
*
*   Determine the maximum occupation of each subshell; override
*   incorrect user input
*
      IF (NORB .GT. NCORE) THEN
         WRITE(istde,*) 'Redefine default maximum relativistic'
     &, ' peel subshell occupations?'
         YES = GETYN ()
      ELSE
         YES = .FALSE.
      ENDIF
*
      IF (YES) THEN
         WRITE(istde,*) 'Enter list of maximum occupations'
     &, ' to be redefined:'
         CALL PRSRCN (NP,SYM,NCORE,NORB,IOCCS)
      ELSE
         DO 4 I = 1,NORB
            IOCCS(I) = 0
    4    CONTINUE
      ENDIF
*
      DO 5 I = 1,NORB
         IF (IOCCS(I) .GT. 0) THEN
            IF (N2J(I) .LE. 7) THEN
               MXOCC(I) = MIN (IOCCS(I),N2J(I)+1)
            ELSE
               MXOCC(I) = 2
            ENDIF
         ELSE
            IF (N2J(I) .LE. 7) THEN
               MXOCC(I) = N2J(I)+1
            ELSE
               MXOCC(I) = 2
            ENDIF
         ENDIF
    5 CONTINUE
*
      IF (YES) THEN
         WRITE(istde,*) 'Maximum peel subshell occupations:'
         CALL PRNTCN (NP,SYM,MXOCC,NCORE,NORB)
      ENDIF
*
*   Determine the minimum occupation of each subshell
*
      IF (NORB .GT. NCORE) THEN
        WRITE(istde,*) 'Redefine default minimum relativistic'
     &, ' peel subshell occupations?'
        YES = GETYN ()
      ELSE
        YES = .FALSE.
      ENDIF
*
      IF (YES) THEN
        WRITE(istde,*) 'Enter list of minimum occupations'
     &, ' to be redefined:'
        CALL PRSRCN (NP,SYM,NCORE,NORB,MNOCC)
      ELSE
        DO 6 I = 1,NORB
          MNOCC(I) = 0
    6   CONTINUE
      ENDIF
*
      IF (YES) THEN
         WRITE(istde,*) 'Minimum peel subshell occupations:'
         CALL PRNTCN (NP,SYM,MNOCC,NCORE,NORB)
      ENDIF
*
*   Read the reference distributions; determine the number of
*   electrons and the atomic parities and J values
*
      NREFD = 2
      CALL ALLOC (PIOREF,NNNW*NREFD,4)
      CALL ALLOC (PNIPAR,       NREFD,4)
      CALL ALLOC (PNTITJ,NNNTJV*NREFD,4)
      CALL ALLOC (PNTNTJ,       NREFD,4)
      NREF = 0
      IF (NORB .GT. NCORE) THEN
    7    WRITE(istde,*) 'Enter a relativistic reference peel'
     &, ' subshell configuration (null if done):'
         CALL PRSRCN (NP,SYM,NCORE,NORB,IOCCS)
         DO 8 I = NCORE+1,NORB
            IF (IOCCS(I) .NE. 0) GOTO 9
    8    CONTINUE
         GOTO 12
    9    LSUM = 0
         ICOUNT = 0
         DO 10 I = NCORE+1,NORB
            IOCCI = IOCCS(I)
            IF ((IOCCI .GT. MXOCC(I)) .OR. (IOCCI .LT. MNOCC(I))) THEN
               WRITE(istde,*) 'GENRL: Subshell occupation out of'
     &, ' bounds; redo ...'
               GOTO 7
            ENDIF
            LSUM = LSUM+IOCCI*NL(I)
            ICOUNT = ICOUNT+IOCCI
            IPTY = (-1)**LSUM
   10    CONTINUE
         IF (ICOUNT .GT. 0) THEN
            IF (NREF .EQ. 0) THEN
               NPEEL = ICOUNT
            ELSE
               IF (ICOUNT .NE. NPEEL) THEN
                  WRITE(istde,*) 'GENRL: Inconsistent number of'
     &, ' peel electrons; redo ...'
                  GOTO 7
               ENDIF
            ENDIF
            NREF = NREF+1
            IF (NREF .GT. NREFD) THEN
               NEWSIZ = NREFD+NREFD/2
               CALL RALC2D (PIOREF,NNNW,NREFD,NNNW,NEWSIZ,4)
               CALL RALLOC (PNIPAR,       NREFD,       NEWSIZ,4)
               CALL RALC2D (PNTITJ,NNNTJV,NREFD,NNNTJV,NEWSIZ,4)
               CALL RALLOC (PNTNTJ,       NREFD,       NEWSIZ,4)
               NREFD = NEWSIZ
            ENDIF
            DO 11 I = NCORE+1,NORB
               IOREF(I,NREF) = IOCCS(I)
   11       CONTINUE
            IPAR(NREF) = IPTY
            IF (IPTY .EQ. -1) THEN
               WRITE(istde,*) 'Parity of configuration is odd'
            ELSE
               WRITE(istde,*) 'Parity of configuration is even'
            ENDIF
            WRITE(istde,*) 'Enter the J values for this configuration'
     &, ' (null for all possible J values):'
            CALL PARSJL (NPEEL,ITJ(1,NREF),NTJ(NREF))
            GOTO 7
         ENDIF
*
   12    CONTINUE
Cww   12    WRITE(istde,*) 'Copy list of relativistic reference'
Cww         WRITE(istde,*) ' configurations to screen?'
Cww         YES = GETYN ()
Cww         IF (YES) THEN
Cww            DO 13 I = 1,NREF
Cww               CALL PRNTCN (NP,SYM,IOREF(1,I),NCORE,NORB)
Cww               CALL PRNTPJ (IPAR(I),NTJ(I),ITJ(1,I))
Cww   13       CONTINUE
Cww         ENDIF
*
*   Determine the number of substitutions allowed
*
   14    WRITE(istde,*) 'Enter the maximum number of substitutions'
         WRITE(istde,*) ' permitted (0, ..., 9):'
         READ (*,'(1I1)',IOSTAT = IOS) NSUB
         IF ((IOS .NE. 0) .OR. (NSUB .LT. 0)) THEN
            WRITE(istde,*) 'GENRL: Unacceptable response;'
     &, ' redo ...'
            GOTO 14
         ENDIF
*
*   Is a copy of the CSF list to appear on screen?
*
Cww         WRITE(istde,*) ' Copy list of CSFs to screen?'
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
      WRITE (21,300) (NP(I),SYM(I),I = 1,NCORE)
      WRITE (21,'(A)') 'Peel subshells:'
      IF ((NORB-NCORE).LE.80) THEN
         WRITE (21,300) (NP(I),SYM(I),I = NCORE+1,NORB)
      ELSE
         WRITE (21,300) (NP(I),SYM(I),I = NCORE+1,NCORE+80)
         WRITE (21,300) (NP(I),SYM(I),I = NCORE+81,NORB)
      ENDIF
      WRITE (21,'(A)') 'CSF(s):'
*
*   Determine the location of the last nonempty subshell in the
*   reference distributions
*
      DO 16 I = NCORE+1,NORB
         DO 15 II = 1,NREF
            IF (IOREF(I,II) .NE. 0) ICHECK = I
   15    CONTINUE
   16 CONTINUE
      ICHECK = ICHECK+1
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
      DO 17 I = 1,NCORE
         IOCCS(I) = MXOCC(I)
   17 CONTINUE
*
      NLEFT = NPEEL
      DO 18 I = NCORE+1,NORB
         IF (NLEFT .GT. 0) THEN
            IOCCS(I) = MIN (MXOCC(I),NLEFT)
            NLEFT = NLEFT-IOCCS(I)
            NLAST = I
         ELSE
            IOCCS(I) = 0
         ENDIF
   18 CONTINUE
*
      IF (NPEEL .GT. 0) THEN
         GOTO 20
      ELSE
         ICOUNT = 1
         GOTO 36
      ENDIF
*
   19 IF (IOCCS(NORB) .EQ. NPEEL) GOTO 36
*
*   Distribute the electrons in the peel subshells
*
      CALL DSTBUT (IOCCS,NCORE,NORB,NLAST)
*
*   Expedite the search by omitting all distributions with
*   more than NSUB electrons outside the last subshell
*   occupied in the reference distributions
*
      IF (IOCCS(ICHECK) .GT. NSUB) THEN
         IOCCS(NORB) = IOCCS(ICHECK)
         IOCCS(ICHECK) = 0
         NLAST = NORB
         GOTO 19
      ENDIF
*
*   Reject this distribution if any peel subshell occupation
*   is out of bounds
*
   20 DO 21 I = NCORE+1,NORB
         IOCT = IOCCS(I)
         IF ((IOCT .GT. MXOCC(I)) .OR.
     :       (IOCT .LT. MNOCC(I))) GOTO 19
   21 CONTINUE
*
*   Compute the parity of this distribution; the parity of the core
*   is always even
*
      LSUM = 0
      DO 22 I = NCORE+1,NORB
         LSUM = LSUM+IOCCS(I)*NL(I)
   22 CONTINUE
      IPTY = (-1)**LSUM
*
*   Reject this distribution if the number of excitations or
*   parity is incorrect with respect to all reference
*   distributions
*
      DO 24 II = 1,NREF
         NEX = 0
         DO 23 I = NCORE+1,NORB
            NEX = NEX + MAX (0,IOREF(I,II)-IOCCS(I))
   23    CONTINUE
         IF (NEX .LE. NSUB) THEN
            IF (IPAR(II) .EQ. IPTY) GOTO 25
         ENDIF
   24 CONTINUE
      GOTO 19
*
*   Determine all peel subshell quantum numbers
*
   25 CALL GETSQN (NCORE,NORB,IOCCS)
*
*   Reset the pointers for the lists of subshell quantum numbers
*
      DO 26 I = NCORE+1,NORB
         IOPTR(I) = 1
   26 CONTINUE
*
*   Determine the first set of all intermediate angular momenta
*
      DO 27 I = 1,NORBM2-NCORE
         ICPTR(I) = 1
         CALL COUPLE (NCORE,I)
   27 CONTINUE
*
*   Initialize level indicator
*
      LVL = NORBM1-NCORE
*
*   Generate all possible final angular momenta
*
   28 DO 29 I = LVL,NORBM1-NCORE
         CALL COUPLE (NCORE,I)
   29 CONTINUE
*
*   Reset level indicator
*
      LVL = NORBM1-NCORE
*
      DO 34 I = 1,LJCL(NORBM1-NCORE)
*
*   Store information only for cases that have the correct J/parity
*   and the correct number of excitations
*
         JATOM = JCLIST(NORBM1-NCORE,I)
         DO 32 II = 1,NREF
            NEX = 0
            DO 30 IV = NCORE+1,NORB
               NEX = NEX + MAX (0,IOREF(IV,II)-IOCCS(IV))
   30       CONTINUE
            IF (NEX .LE. NSUB) THEN
               IF (IPAR(II) .EQ. IPTY) THEN
                  IF (NTJ(II) .EQ. -1) GOTO 33
                  DO 31 III = 1,NTJ(II)
                     IF (ITJ(III,II) .EQ. JATOM) GOTO 33
   31             CONTINUE
               ENDIF
            ENDIF
   32    CONTINUE
*
*   Match not found; skip to end of loop
*
         GOTO 34
*
*   Match found
*
   33    ICOUNT = ICOUNT+1
         CALL WRTCSF (NCORE,NORB,IOCCS,JATOM,IPTY,PRNTSC)
*
   34 CONTINUE
*
*   Reset the pointers and determine the level indicator
*
   35 LVLP1 = LVL+1
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
               GOTO 35
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
           IF (IOPTR(NCORE+1) .GT. LJL(NCORE+1)) GOTO 19
*
         ENDIF
      ENDIF
*
*   Obtain a new set of couplings with the same distribution
*
      GOTO 28
*
*   Print summary
*
   36 WRITE(istde,*) NORB,' relativistic subshells.'
     &, ICOUNT,' relativistic CSFs generated.'
*
*   All done
*
      CALL DALLOC (PIOREF)
      CALL DALLOC (PNIPAR)
      CALL DALLOC (PNTITJ)
      CALL DALLOC (PNTNTJ)
      CALL DALLOC (PJCLST)
*
      RETURN
*
  300 FORMAT (120(1X,1I2,1A2))
*
      END
