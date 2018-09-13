************************************************************************
*                                                                      *
      SUBROUTINE DENSMCP(DINT1,DINT2,DINT3,DINT4,DINT5,DINT6)       
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
      DIMENSION DINT1(NNNW,NNNW) 
      DIMENSION DINT2(NNNW,NNNW) 
      DIMENSION DINT3(NNNW,NNNW) 
      DIMENSION DINT4(NNNW,NNNW) 
      DIMENSION DINT5(NNNW,NNNW) 
      DIMENSION DINT6(NNNW,NNNW) 
*
      POINTER (PCOEFF,COEFF(*))
      POINTER (PICLMN,ICLMN(*))
      POINTER (PINDEX,INDEX(*))
*
      POINTER (PNTEMT1,EMT1(*))
      POINTER (PNTEMT2,EMT2(*))
      POINTER (PNTEMT3,EMT3(*))
      POINTER (PNTEMT4,EMT4(*))
      POINTER (PNTEMT5,EMT5(*))
      POINTER (PNTEMT6,EMT6(*))
      POINTER (PIENDC,IENDC(0:*))
      POINTER (PNIROW,IROW(*))
*
      POINTER (PNEVEC,EVEC(*))
*
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
*   Allocate storage that is local to this subroutine
*
      NDIM = 1
      CALL ALLOC (PCOEFF,NDIM,8)
      CALL ALLOC (PICLMN,NDIM,4)
      CALL ALLOC (PINDEX,NDIM,4)
*
      READ (30) NELMNT                                                  
      CALL ALLOC (PIENDC,NCF+1,4)                                       
      CALL ALLOC (PNIROW,NELMNT,4)
      READ (30) (IENDC(I),I = 0,NCF),(IROW(I),I = 1,NELMNT)
      CLOSE (30)
*
*   Other initializations
*
      CALL ALLOC (PNTEMT1,NELMNT,8)
      CALL ALLOC (PNTEMT2,NELMNT,8)
      CALL ALLOC (PNTEMT3,NELMNT,8)
      CALL ALLOC (PNTEMT4,NELMNT,8)
      CALL ALLOC (PNTEMT5,NELMNT,8)
      CALL ALLOC (PNTEMT6,NELMNT,8)
*
      NWM1 = NW-1
*
      DO 1 I = 1,NELMNT
         EMT1(I) = 0.0D 00
         EMT2(I) = 0.0D 00
         EMT3(I) = 0.0D 00
         EMT4(I) = 0.0D 00
         EMT5(I) = 0.0D 00
         EMT6(I) = 0.0D 00
    1 CONTINUE
*   
*   Accumulate diagonal terms that do not require MCP coefficients
*   
*   
*   Piece involving I(a,a) integrals
*   
      DO 3 IA = 1,NW
         SET = .FALSE.
         DO 2 IR = 1,NCF
            QA = DBLE (IQ (IA,IR))
            IF (QA .GT. 0.0D 00) THEN
               IF (.NOT. SET) THEN
                  DIAA1 = DINT1(IA,IA)
                  DIAA2 = DINT2(IA,IA)
                  DIAA3 = DINT3(IA,IA)
                  DIAA4 = DINT4(IA,IA)
                  DIAA5 = DINT5(IA,IA)
                  DIAA6 = DINT6(IA,IA)
                  SET = .TRUE. 
               ENDIF
               IDIAG = IENDC(IR-1)+1
               EMT1(IDIAG) = EMT1(IDIAG)+QA*DIAA1
               EMT2(IDIAG) = EMT2(IDIAG)+QA*DIAA2
               EMT3(IDIAG) = EMT3(IDIAG)+QA*DIAA3
               EMT4(IDIAG) = EMT4(IDIAG)+QA*DIAA4
               EMT5(IDIAG) = EMT5(IDIAG)+QA*DIAA5
               EMT6(IDIAG) = EMT6(IDIAG)+QA*DIAA6
            ENDIF
    2    CONTINUE
    3 CONTINUE
*                     
*   Accumulate one-electron terms that require MCP coefficients
*                                    
      REWIND (31)                         
      READ (31)                                           
      READ (31)                   
      READ (31)                             
*                                                     
*   Attempt to read another block of data  
*                                                       
   16 READ (31,IOSTAT = IOS) LAB,NCONTR      
*                                 
      IF (IOS .EQ. 0) THEN                     
*                                                        
*   Read successful; decode the labels of I(ab) 
*
         IA = MOD (LAB,KEY)
         IB = LAB/KEY
*  
*   Compute I(ab)
*  
         TEGRAL1 = DINT1(IA,IB)
         TEGRAL2 = DINT2(IA,IB)
         TEGRAL3 = DINT3(IA,IB)
         TEGRAL4 = DINT4(IA,IB)
         TEGRAL5 = DINT5(IA,IB)
         TEGRAL6 = DINT6(IA,IB)
*  
*   Ensure that storage is adequate to read in the rest of
*   this block
*  
         IF (NCONTR .GT. NDIM) THEN
            CALL DALLOC (PCOEFF)
            CALL DALLOC (PICLMN)
            CALL DALLOC (PINDEX)
            NDIM = NCONTR
            CALL ALLOC (PCOEFF,NDIM,8)
            CALL ALLOC (PICLMN,NDIM,4)
            CALL ALLOC (PINDEX,NDIM,4)
         ENDIF
*   
*   Read the column index, the sparse matrix index, and the
*   coefficient for all contributions from this integral
*
         READ (31) (ICLMN(I),INDEX(I),COEFF(I),I = 1,NCONTR)
*  
*   Store all the contributions from this integral
*  
         DO 17 I = 1,NCONTR
            LOC = INDEX(I)
            EMT1(LOC) = EMT1(LOC)+TEGRAL1*COEFF(I)
            EMT2(LOC) = EMT2(LOC)+TEGRAL2*COEFF(I)
            EMT3(LOC) = EMT3(LOC)+TEGRAL3*COEFF(I)
            EMT4(LOC) = EMT4(LOC)+TEGRAL4*COEFF(I)
            EMT5(LOC) = EMT5(LOC)+TEGRAL5*COEFF(I)
            EMT6(LOC) = EMT6(LOC)+TEGRAL6*COEFF(I)
   17    CONTINUE
*  
*   Return to the start of the loop
*   
         GOTO 16
*   
      ENDIF
*
*   Deallocate storage that is local to this routine
*
      CALL DALLOC (PCOEFF)
      CALL DALLOC (PICLMN)
      CALL DALLOC (PINDEX)

      ICI = 0
      DO 21 I = 1,NELMNT
         IRI = IROW(I)
         IF (I .GT. IENDC(ICI)) ICI = ICI+1
         DO 22 J = 1,NVEC
            LOC = (J-1)*NCF
            CONTRI1 = EVEC(ICI+LOC)*EVEC(IRI+LOC)*EMT1(I)
            CONTRI2 = EVEC(ICI+LOC)*EVEC(IRI+LOC)*EMT2(I)
            CONTRI3 = EVEC(ICI+LOC)*EVEC(IRI+LOC)*EMT3(I)
            CONTRI4 = EVEC(ICI+LOC)*EVEC(IRI+LOC)*EMT4(I)
            CONTRI5 = EVEC(ICI+LOC)*EVEC(IRI+LOC)*EMT5(I)
            CONTRI6 = EVEC(ICI+LOC)*EVEC(IRI+LOC)*EMT6(I)
            IF (IRI.NE.ICI) THEN
               CONTRI1 = 2.0D 00 * CONTRI1
               CONTRI2 = 2.0D 00 * CONTRI2
               CONTRI3 = 2.0D 00 * CONTRI3
               CONTRI4 = 2.0D 00 * CONTRI4
               CONTRI5 = 2.0D 00 * CONTRI5
               CONTRI6 = 2.0D 00 * CONTRI6
            ENDIF  
            DENS1(J) = DENS1(J) + CONTRI1
            DENS2(J) = DENS2(J) + CONTRI2
            DENS3(J) = DENS3(J) + CONTRI3
            DENS4(J) = DENS4(J) + CONTRI4
            DENS5(J) = DENS5(J) + CONTRI5
            DENS6(J) = DENS6(J) + CONTRI6
   22    CONTINUE 
   21 CONTINUE
      CALL DALLOC (PNTEMT1)
      CALL DALLOC (PNTEMT2)
      CALL DALLOC (PNTEMT3)
      CALL DALLOC (PNTEMT4)
      CALL DALLOC (PNTEMT5)
      CALL DALLOC (PNTEMT6)
      CALL DALLOC (PIENDC)
      CALL DALLOC (PNIROW)
*
*   Close the angular files
*
      DO 9 I = 30,32+KMAXF
        CLOSE (I)
    9 CONTINUE

      RETURN
      END
