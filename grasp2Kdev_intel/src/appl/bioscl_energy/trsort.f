************************************************************************
*                                                                      *
      SUBROUTINE TRSORT (NAME,NFILE,NFILE2,LPRINT,JKP,IBLKI,IBLKF)
*                                                                      *
*   Routine to sort angular coefficients into list based on integral   *
*   labels rather than CSF.  A tree sort is used. To save space, the   *
*   CSF pair labels and coefficients  are not read into arrays until   *
*   the sorting has been done.                                         *
*                                                                      *
*   Call(s) to: [OSCL92]: ALCLLA, ALCNMA.                              *
*               [LIB92]: ALLOC, DALLOC.                                *
*                                                                      *
*   Original author(s) unknown. Modifications for dynamic storage by   *
*   Farid A. Parpia.                                                   *
*                                                                      *
*                                           Last update: 28 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)

      include 'parameters.def'
CFF      PARAMETER (NNNW = 120)
CFF      PARAMETER (NCA = 65536,
CFF     :           KEY = 121)
      PARAMETER (KEY = KEYORB)
      PARAMETER (NCA = 65536)

Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL FIRST,LPRINT
      CHARACTER*2 S(-9:9)
      CHARACTER*2 NH
      CHARACTER*24 NAME(2)
*
      POINTER (PIBEG,IBEG(1)),(PILAB,ILAB(1)),(PILAST,ILAST(1)),
     :        (PILEFT,ILEFT(1)),(PIPTCS,IPTCSF(1)),(PIRIGH,IRIGHT(1)),
     :        (PLBLIN,LBLINT(1))
      POINTER (PIPTR,IPTR(1)),(PISLDR,ISLDR(1)),(PXSLDR,XSLDR(1))
      POINTER (PISLDR1,ISLDR1(1))
*
      POINTER (PLABEL,JLABL(1)),(PCOEFF,XL(1))
      POINTER (PNTRKP,KP(1))
*
      COMMON/ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /OSC6/NKP,PNTRKP
     :      /DEFAULT/NDEF,NDUMP
*
      S(-9)  = '-9'
      S(-8)  = '-8'                                                     
      S(-7)  = '-7'                                                     
      S(-6)  = '-6'                                                     
      S(-5)  = '-5'                                                     
      S(-4)  = '-4'                                                     
      S(-3)  = '-3'                                                     
      S(-2)  = '-2'                                                     
      S(-1)  = '-1'                                                     
      S(0)   = '+0'                                                     
      S(1)   = '+1'   
      S(2)   = '+2'                
      S(3)   = '+3'       
      S(4)   = '+4'   
      S(5)   = '+5'       
      S(6)   = '+6'    
      S(7)   = '+7'
      S(8)   = '+8'
      S(9)   = '+9'
*
*   Position file at beginning of list of integrals
*
      REWIND (NFILE)
*
*   Initialize
*
      FIRST = .TRUE.
      NMCP = 0
      NINT = 0
*
*   Initial allocation of storage to local arrays
*
      NLABEL = 1
      CALL ALLOC (PLABEL,NLABEL,4)
      CALL ALLOC (PCOEFF,NLABEL,8)
*
      CALL ALCLLA (PIBEG,PILAB,PILAST,PILEFT,
     :             PIPTCS,PIRIGH,PLBLIN,LLDIM,1)
      CALL ALCNMA (PIPTR,PISLDR,PISLDR1,PXSLDR,NMDIM,1)
*
*   Now the rest of the elements
*
    1 READ (NFILE,END = 12) IR,IS,NI
      IF (NI .GT. NLABEL) THEN
         CALL RALLOC (PLABEL,NLABEL,NI,4)
         CALL RALLOC (PCOEFF,NLABEL,NI,8)
         NLABEL = NI
      ENDIF
      READ (NFILE,END = 99,ERR = 99) (JLABL(I),XL(I),I = 1,NI)
      IF ((IR .EQ. 0) .OR. (IS .EQ. 0) .OR. (NI .EQ. 0)) GOTO 1
*
      IF (FIRST) THEN
*
*   List is empty
*
         M = 0
         J = 0
*
*   Set up list pointers and insert first element
*
         ICOUNT = 0
    3    ICOUNT = ICOUNT+1
         IF (ICOUNT .GT. NI) GOTO 1
         IF (JLABL(ICOUNT) .EQ. 0) GOTO 3
         ILAB(1) = JLABL(ICOUNT)
*
         IRIGHT(1) = 0
         ILEFT(1) = 0
         IBEG(1) = 1
         IPTR(1) = 0
         ILAST(1) = 0
*
         M = 1
         J = 1
*
         FIRST = .FALSE.
*
      ELSE
*
         ICOUNT = 0
*
      ENDIF
*
*   Sort integral list using tree sort
*
*   Take next nonzero element
*
    4 ICOUNT = ICOUNT+1
      IF (ICOUNT .GT. NI) GOTO 1
      JLAB = JLABL(ICOUNT)
      IF (JLAB .EQ. 0) GOTO 4
      X = XL (ICOUNT)
*
      M = M+1
      IF (M .GT. NMDIM) CALL ALCNMA(PIPTR,PISLDR,PISLDR1,PXSLDR,NMDIM,2)
      I = 1
*
*   Search for place in tree
*
    5 IF (JLAB-ILAB(I)) 6,10,8
    6 K = IRIGHT(I)
      IF (K .NE. 0) GOTO 7
      J = J+1
      IF (J .GT. LLDIM) CALL ALCLLA (PIBEG,PILAB,PILAST,PILEFT,
     :                               PIPTCS,PIRIGH,PLBLIN,LLDIM,2)
      IRIGHT(I) = J
      GOTO 9
    7 I = K
      GOTO 5
    8 K = ILEFT(I)
      IF (K .NE. 0) GOTO 7
      J = J+1
      IF (J .GT. LLDIM) CALL ALCLLA (PIBEG,PILAB,PILAST,PILEFT,
     :                               PIPTCS,PIRIGH,PLBLIN,LLDIM,2)
      ILEFT(I) = J
*
*   When found, update list.
*
    9 ILAST(J) = I
      IRIGHT(J) = 0
      ILEFT(J) = 0
      IBEG(J) = M
      ILAB(J) = JLAB
      IPTR(M) = 0
      GOTO 4
   10 K = IBEG(I)
   11 L = K
      K = IPTR(L)
      IF (K .NE. 0) GOTO 11
      IPTR(L) = M
      IPTR(M) = 0
      GOTO 4
*
*   The end of the CSF-based file has been reached
*
   12 IF (FIRST .OR. (M .EQ. 0)) GOTO 20
*
*  Sort is complete. Unpack list
*
      NMCP = M
      NINT = J
      L = 0
      M = 0
      I = 1
*
*   Search for smallest element
*
   13 K = IRIGHT(I)
      IF (K .EQ. 0) GOTO 14
      I = K
      GOTO 13
*
*   Insert in sorted list
*
   14 IF (ILAB(I) .EQ. 0) GOTO 16
      L = L+1
      LBLINT(L) = ILAB(I)
      K = IBEG(I)
*
*   Copy list of pointers to CSF/coefficients into new list
*
   15 M = M+1
      ISLDR(M) = K
      K = IPTR(K)
      IF (K .NE. 0) GOTO 15
      IPTCSF(L) = M
      ILAB(I) = 0
*
*   Next smallest element is on left of last element
*
      K = ILEFT(I)
      IF (K .EQ. 0) GOTO 16
      I = K
      GOTO 13
*
*   If no element on left, next smallest is previous element
*
   16 I = ILAST(I)
      IF (I .NE. 0) GOTO 14
*
*   List is unpacked. Invert CSF/coefficient pointer list to give
*   position list for CSF/coefficients as they are read in
*
      DO 17 I = 1,NMCP
         K = ISLDR(I)
         IPTR(K) = I
   17 CONTINUE
*
*   Now read CSF pairs and coefficients into correct positions in
*   sorted list
*
      IMCP = 0
      REWIND (NFILE)
   18 READ (NFILE,END = 20) IR,IS,NI
      IF ((IR .EQ. 0) .OR. (IS .EQ. 0) .OR. (NI .EQ. 0)) GOTO 18
      READ (NFILE) (JLABL(I),XL(I), I = 1,NI)
      DO 19 I = 1,NI
         IF (JLABL(I) .NE. 0) THEN
            IMCP = IMCP+1
            K = IPTR(IMCP)
            ISLDR(K) = IR
            ISLDR1(K) = IS
c           ISLDR(K) = IR*NCA+IS
            XSLDR(K) = XL(I)
         ENDIF
   19 CONTINUE
      GOTO 18
*
*   The integral-based list is completely known
*
   20 REWIND (NFILE)
*                                                                       
*  If first set of data open the file and print                         
*  some data to later be able to identify the file                      
*
      if(iblki.eq.1.and.iblkf.eq.1) then
        J1 = INDEX(NAME(1),' ')
        J2 = INDEX(NAME(2),' ')
        OPEN (UNIT = NFILE2,
     :    FILE=NAME(1)(1:J1-1)//'.'//NAME(2)(1:J2-1)//'.'
     :    //S(KP(JKP))//'T',
     :          STATUS='UNKNOWN',FORM='UNFORMATTED')                    
      endif
        WRITE (NFILE2) IBLKI,IBLKF,NW,NKP
         WRITE (NFILE2) NINT
      IF (NMCP .NE. 0) THEN
         IF (LPRINT) WRITE (99,301)
         MLR = 1
         DO 23 I = 1,NINT
            MUP = IPTCSF(I)
            INTS = MUP-MLR+1
               WRITE (NFILE2) LBLINT(I),INTS
            IF (LPRINT) THEN
               IA = LBLINT(I)/KEY
               IB = MOD (LBLINT(I),KEY)
               WRITE (99,302) NP(IA),NH(IA),NP(IB),NH(IB)
            ENDIF
               WRITE (NFILE2) (ISLDR(M),ISLDR1(M),XSLDR(M),M = MLR,MUP)
            IF (LPRINT) THEN
               DO 22 M = MLR,MUP
c                  IS = MOD (ISLDR(M),NCA)
                  IS = ISLDR1(M)
c                 IR = ISLDR(M)/NCA
                  IR = ISLDR(M)
                  WRITE (99,303) IR,IS,XSLDR(M)
   22          CONTINUE
            ENDIF
            MLR = MUP+1
   23    CONTINUE
*
      ENDIF
*
*   Deallocate storage for local arrays
*
      CALL DALLOC (PLABEL)
      CALL DALLOC (PCOEFF)
      CALL ALCLLA (PIBEG,PILAB,PILAST,PILEFT,
     :             PIPTCS,PIRIGH,PLBLIN,LLDIM,3)
      CALL ALCNMA (PIPTR,PISLDR,PISLDR1,PXSLDR,NMDIM,3)
*
      RETURN
*
*   Error handling
*
   99 PRINT *, 'TRSORT: Error reading CSF-based file.'
      STOP
*
  301 FORMAT (//'  k'
     :         /' d  (rs)  Coefficients:'
     :         /'  ab'
     :       ///'    a     b          r         s  Coefficient'/)
  302 FORMAT (2(2X,I2,A2))
  303 FORMAT (14X,1I6,2X,1I6,2X,1P,1D22.15)
*
      END
