************************************************************************
*                                                                      *
      SUBROUTINE SMSMCP(VINT)       
*                                                                      *
*   This routine controls the main sequence of routine calls for the   *
*   calculation  of the  sms parameter, the electron density at the    *
*   origin and the field shift between two isotopes characterized      *
*   by fermi distributions c and a                                     *
*                                                                      *
*   Call(s) to: [LIB92]: ALCBUF, ALLOC, CONVRT, DALLOC, GETYN          *
*                        ITJPO, NUCPOT, RKCO, TNSRJJ,                  *
*               [SMS92]: RINTISO, RINTDENS, VINTI                      *
*                                                                      *
*   Written by Per Jonsson                                             *
*                                                                      *
*                                         Last revision: 10 Nov 1995   *
*                                                                      *
************************************************************************
*     
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PNEVAL,PNTRIQ,PNIVEC
      POINTER (PNEVAL,EVALDUMMY)
      POINTER (PNTRIQ,RIQDUMMY)
      POINTER (PNIVEC,IVECDUMMY)

      LOGICAL SET
      DIMENSION VINT(NNNW,NNNW) 
*
      POINTER (PCOEFF,COEFF(*))
      POINTER (PICLMN,ICLMN(*))
      POINTER (PINDEX,INDEX(*))
      POINTER (PNSWAP,NSWAP(*))
*
      POINTER (PNTEMT,EMT(*))
      POINTER (PIENDC,IENDC(0:*))
      POINTER (PNIROW,IROW(*))
*
      POINTER (PNEVEC,EVEC(*))
*
cff      POINTER (PNDENS,DENS(*))
      POINTER (PNSMS,SMSC(*))
      POINTER (PNDENS1,DENS1(*))
      POINTER (PNDENS2,DENS2(*))
      POINTER (PNDENS3,DENS3(*))
      POINTER (PNDENS4,DENS4(*))
      POINTER (PNDENS5,DENS5(*))
      POINTER (PNDENS6,DENS6(*))
*
      COMMON/EIGVAL/EAV,PNEVAL
     :      /EIGVEC/PNEVEC
     :      /HMAT/PNTEMT,PIENDC,PNIROW,NELMNT
     :      /MCPA/KMAXF
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
     :      /PRNT/NVEC,PNIVEC,NVECMX
     :      /SMS1/PNSMS,PNDENS1,PNDENS2,PNDENS3,PNDENS4,PNDENS5,PNDENS6
*
CGG      PARAMETER (KEYORB = 121)
      PARAMETER (KEY = KEYORB)
*
      OPEN (UNIT=20,FILE='sms.20',FORM='UNFORMATTED',                   
     :      STATUS='OLD')     

*
*   Allocate storage that is local to this subroutine
*
      NDIM = 1
      CALL ALLOC (PCOEFF,NDIM,8)
      CALL ALLOC (PICLMN,NDIM,4)
      CALL ALLOC (PINDEX,NDIM,4)
      CALL ALLOC (PNSWAP,NDIM,4)
*
*   General initializations
*
      READ (30) NELMNT
      CALL ALLOC (PIENDC,NCF+1,4)
      CALL ALLOC (PNIROW,NELMNT,4)
      READ (30) (IENDC(I),I = 0,NCF),(IROW(I),I = 1,NELMNT)
      CLOSE (30)
*
*   Other initializations
*
      CALL ALLOC (PNTEMT,NELMNT,8)
*
      NWM1 = NW-1
*
      DO 1 I = 1,NELMNT
         EMT(I) = 0.0D 00
    1 CONTINUE
*   
*   Accumulate diagonal terms that do not require MCP coefficients
*   
*                    k
*   Piece involving G (a,b) integrals
*
      KM = 0
      K = 1
      DO 15 IA = 1,NWM1
         NKJIA = NKJ(IA)
         IAP1 = IA+1
         DO 14 IB = IAP1,NW
            NKJIB = NKJ(IB)
            SET = .FALSE.
            IF (NAK(IA)*NAK(IB) .GT. 0) THEN
               KMIN = ABS ((NKJIA-NKJIB)/2)
            ELSE
               KMIN = ABS ((NKJIA-NKJIB)/2)+1
            ENDIF
            IF (MOD (K-KMIN,2) .EQ. 0) THEN
               KMAX = (NKJIA+NKJIB)/2
               IF (KMAX .GT. KM) KM = KMAX
               IF ((K .GE. KMIN) .AND. (K .LE. KMAX)) THEN
                  DO 13 IR = 1,NCF
                     COEF = GCO (K,IR,IA,IB)
                     IF (ABS (COEF) .GT. 0.0D 00) THEN
                        IF (.NOT. SET) THEN
Cww
                           GKAB = VINT (IA,IB)*VINT (IB,IA)
Cww
                           SET = .TRUE.
                        ENDIF
                        IDIAG = IENDC(IR-1)+1
                        EMT(IDIAG) = EMT(IDIAG)-COEF*GKAB
                     ENDIF
   13             CONTINUE
               ENDIF
            ENDIF
   14    CONTINUE
   15 CONTINUE
*
*   Accumulate two-electron terms that require MCP coefficients
*
      NFILE = 33
*
      REWIND (NFILE)
      REWIND (20)
      READ (NFILE)
      READ (NFILE)
      READ (NFILE)
*
*   The multipolarity of the integral can be deduced from the file
*   unit number
*
      K = NFILE-32
*
*   Attempt to read another block of data
*
   18 READ (NFILE,IOSTAT = IOS) LAB,NCONTR
*
         IF (IOS .EQ. 0) THEN
*                                          k
*   Read successful; decode the labels of R (abcd)
*
            ID = MOD (LAB,KEY)
            LAB = LAB/KEY
            IB = MOD (LAB,KEY)
            LAB = LAB/KEY
            IC = MOD (LAB,KEY)
            IA = LAB/KEY
*
*   Compute the Vinti integrals
*
Cww
Cww            TEGRAL = VINT (IA,IB)*VINT (IC,ID)
            TEGRAL = VINT (IA,IC)*VINT (IB,ID)
Cww
*
*   Ensure that storage is adequate to read in the rest of
*   this block
*
            IF (NCONTR .GT. NDIM) THEN
               CALL DALLOC (PCOEFF)
               CALL DALLOC (PICLMN)
               CALL DALLOC (PINDEX)
               CALL DALLOC (PNSWAP)
               NDIM = NCONTR
               CALL ALLOC (PCOEFF,NDIM,8)
               CALL ALLOC (PICLMN,NDIM,4)
               CALL ALLOC (PINDEX,NDIM,4)
               CALL ALLOC (PNSWAP,NDIM,4) 
            ENDIF
*
*   Read the column index, the sparse matrix index, and the
*   coefficient for all contributions from this integral
*
            READ (NFILE) (ICLMN(I),INDEX(I),COEFF(I),I = 1,NCONTR)
            READ (20) (NSWAP(I),I = 1,NCONTR)
*
*   Store all the contributions from this integral
*
            DO 19 I = 1,NCONTR
               LOC = INDEX(I)
               EMT(LOC) = EMT(LOC)-TEGRAL*COEFF(I)*(-1)**NSWAP(I)
   19       CONTINUE
*
*   Return to the start of the loop
*
            GOTO 18
*
         ENDIF
*
*
*   Deallocate storage that is local to this routine
*
      CALL DALLOC (PCOEFF)
      CALL DALLOC (PICLMN)
      CALL DALLOC (PINDEX)
      CALL DALLOC (PNSWAP)

      ICI = 0
      DO 21 I = 1,NELMNT
         IRI = IROW(I)
         IF (I .GT. IENDC(ICI)) ICI = ICI+1
         DO 22 J = 1,NVEC
            LOC = (J-1)*NCF
            CONTRI = EVEC(ICI+LOC)*EVEC(IRI+LOC)*EMT(I)
            IF (IRI.NE.ICI) THEN
               CONTRI = 2.0D 00 * CONTRI
            ENDIF  
            SMSC(J) = SMSC(J) + CONTRI
   22    CONTINUE 
   21 CONTINUE
      CALL DALLOC (PNTEMT)
      CALL DALLOC (PIENDC)                                                   
      CALL DALLOC (PNIROW)
*
*   Close the angular files
*
      CLOSE (20)
      DO 89 I = 30,32+KMAXF
        CLOSE (I)
   89 CONTINUE
      RETURN
      END
