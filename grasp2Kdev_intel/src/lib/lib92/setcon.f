************************************************************************
*                                                                      *
      SUBROUTINE SETCON
*                                                                      *
*   This  subprogram  sets the values of the fundamental and derived   *
*   physical constants, and other useful constants.                    *
*                                                                      *
*   Written by Farid A Parpia, at Oxford   Last updated: 06 Oct 1992   *
*   The values of the physical constants   Last updated: 22 July 2012  *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      LOGICAL LDBPG
*
      COMMON/DEBUGG/LDBPG(5)
     :      /DEF3/EMPAM,RBCM
     :      /DEF9/CVAC,PI
     :      /DEF10/AUCM,AUEV,CCMS,FASI,FBSI
     :      /DEF11/FMTOAU,AUMAMU
*
*   Physical constants from http://physics.nist.gov (2012)
*   Physical constants: see description below.
*
      PARAMETER (AINFCM =       0.529 177 210 92  D-08,
     :           ALFAI  =     137.035 999 074  D 00,
     :           CCMPS  =       2.99 792 458   D 10,
     :           EESU   =       4.803 204 506  D-10,
     :           EMEG   =       9.109 382 91   D-28,
     :           EMEAMU =       5.485 799 0946 D-04,
     :           EMPAMU =       1.007 276 466 812  D 00,
     :           HBARES =       1.054 571 726  D-27,
     :           RINFEV =      13.605 692 53   D 00,
     :           RINFK  = 109 737.315 685 39   D 00)
*
*   Previous values of physical constants are from:
*    (see setcon.f-????)
*   The values of the physical constants used here are taken from
*   E R Cohen and B N Taylor, The 1986 Adjstment of the Fundamental
*   Physical Constants, Report of the CODATA Task Group on Funda-
*   mental COnstants, CODATA Bulletin 63, Pergamon, Elmsford, NY
*   (1986)
*
*      AINFCM: Bohr radius in cm
*      ALFAI : Inverse of the fine-structure constant
*      CCMPS : Speed of light in cm/s
*      EESU  : Electron charge in esu
c              conversion : EESU = EEC * 2997924580 
*      EMEG  : Electron mass in g
*      EMEAMU: Electron mass in amu
*      EMPAMU: Proton mass in amu
*      HBARES: Rationalized Planck's constant in erg s
*      RINFEV: Rydberg in eV
*      RINFK : Rydberg in Kaysers
*
*   In the context of new results, Taylor warns that ... since the
*   output values of a least-squares adjustment are related in a
*   complex way and a change in the measured value of one constant
*   usually leads to corresponding changes in the adjusted values
*   of others, one must be cautious in carrying out calculations
*   using both the [above values] and the results of more recent
*   measurements.
*
*   Calculate constants for /DEF3/:
*
      EMPAM = EMPAMU
      RBCM = AINFCM
*
*   Calculate constants for  /DEF11/:
*
*      CVAC is the speed of light in the vacuum in atomic units
*      B1 is the conversion factor between the amu and the atomic
*         unit of mass
*
      CVAC = ALFAI
      AUMAMU = EMEAMU
*
*   Calculate constants for  /DEF10/:
*
*      AUCM converts au to kaysers assuming an infinitely-heavy
*           nucleus
*      AUEV Converst au to eV assuming an infinitely-heavy
*           nucleus
*      CCMS is the speed of light in the vacuum in centimetres
*            per second
*      FASI converts the Einstein A coefficients from atomic to SI
*            units
*      FBSI converts the Einstein B coefficients from atomic to SI
*            units
*
      AUCM = 2.0D 00*RINFK
      AUEV = 2.0D 00*RINFEV
      CCMS = CCMPS
      FASI = (EMEG/HBARES)*((EESU*EESU/HBARES)**2)
      FBSI = (10.0D 00*AINFCM**3/HBARES)*FASI
*
*   Calculate conversion factor for Fermis to Bohr radii
*
      FMTOAU = 1.0D-13/AINFCM
*
*   Calculate \pi from FORTRAN function
*
      PI = 4.0D 00*ATAN (1.0D 00)
*
*   Printouts
*
      IF (LDBPG(2)) WRITE (99,300) AINFCM,ALFAI,CCMPS,EESU,
     :                             EMEG,EMEAMU,EMPAMU,HBARES,
     :                             RINFEV,RINFK
*
      RETURN
*
  300 FORMAT (/'From SUBROUTINE SETCON:'
     :       /' AINFCM (Bohr radius in cm): ',0P,1D15.9,','
     :       /' ALFAI (Inverse of the fine-structure constant): ',
     :             3P,1D15.9,','
     :       /' CCMPS (Speed of light in cm/s): ',1P,1D14.8,','
     :       /' EESU (Electron charge in esu): ' ,1P,1D14.8,','
     :       /' EMEG (Electron mass in g): ',1P,1D13.7,','
     :       /' EMEAMU (Electron mass in u): ',1P,1D15.9,','
     :       /' EMPAMU (Proton mass in u): ',1P,1D15.9,','
     :       /' HBARES (Rationalized Planck constant in erg s): ',
     :             1P,1D14.8,','
     :       /' RINFEV (Rydberg in eV): ',2P,1D14.8,','
     :       /' RINFK (Rydberg in Kaysers): ',6P,1D16.10,'.')
*
      END
