************************************************************************
*                                                                      *
      BLOCK DATA CONSTS
*                                                                      *
*   Sets up commonly used numerical constants.                         *
*                                                                      *
*                                          Last updated: 06 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
*
      COMMON/CONS/ZERO,HALF,TENTH,ONE,TWO,THREE,TEN
*
      DATA ZERO,HALF,TENTH,ONE,TWO,THREE,TEN/
     :     0.0D 00,0.5D 00,0.1D 00,
     :     1.0D 00,2.0D 00,3.0D 00,
     :    10.0D 00/
*

! Standard input, output and error unit numbers
! CRAY also uses 0 for error, most others use 0
!XHH
      COMMON/iounit/istdi,istdo,istde
      DATA istdi,istdo,istde/5,6,0/

      END
************************************************************************
*                                                                      *
      BLOCK DATA JLABEL
*                                                                      *
*   Sets up strings giving J-values and parity signs. The string for   *
*   angular momentum J is that in the array location 2J+1. If PARITY   *
*   is stored as -1 for odd parity and 1 for even parity, the appro-   *
*   priate array element is in location (3+PARITY)/2.                  *
*                                                                      *
*                                          Last update: 06 Oct 1992    *
*                                                                      *
************************************************************************
*
      CHARACTER*4 JLBL,JLBR,JLBP
*
      COMMON/JLABL/JLBL(32),JLBR(32),JLBP(2)
*
*   Left-justified strings
*
      DATA JLBL/'0   ','1/2 ','1   ','3/2 ','2   ','5/2 ','3   ','7/2 ',
     :          '4   ','9/2 ','5   ','11/2','6   ','13/2','7   ','15/2',
     :          '8   ','17/2','9   ','19/2','10  ','21/2','11  ','23/2',
     :          '12  ','25/2','13  ','27/2','14  ','29/2','15  ','31/2'/
*
*   Right-justified strings
*
      DATA JLBR/'   0',' 1/2','   1',' 3/2','   2',' 5/2','   3',' 7/2',
     :          '   4',' 9/2','   5','11/2','   6','13/2','   7','15/2',
     :          '   8','17/2','   9','19/2','  10','21/2','  11','23/2',
     :          '  12','25/2','  13','27/2','  14','29/2','  15','31/2'/
*
*   Parity signs
*
      DATA JLBP/' -  ',' +  '/
*
      END
************************************************************************
*                                                                      *
      BLOCK DATA TERM
*                                                                      *
*   Taken from GRASP2.                                                 *
*                                                                      *
*   Table of subshell quantum numbers.    Symmetry  of the table for   *
*   particle/hole configurations is used to compress it.               *
*                                                                      *
*                                          Last updated: 30 Sep 1992   *
*                                                                      *
************************************************************************
*
      COMMON/TERMS/NROWS,ITAB(31),JTAB(32),NTAB(327)
*
      DATA NROWS/31/
*
*   A row is defined by a subshell angular momentum and an occupation
*   number
*
*   Each entry ITAB gives the number of terms in a row
*
*   Each entry JTAB gives the starting location -1 of the first triad
*   in a row
*
*   Each triad in NTAB is (v,w,2J+1); here v is the seniority,
*   w resolves any degeneracy in the seniority scheme, and J is the
*   subshell total angular momentum
*
*   Empty subshell or full subshell
*
      DATA (ITAB(I),I =   1,  1)/
     :  1/
      DATA (JTAB(I),I =   1,  1)/
     :  0/
      DATA (NTAB(I),I =   1,  3)/
     :  0,0, 1/
*
*   s, p-   (j = 1/2)
*
      DATA (ITAB(I),I =   2,  2)/
     :  1/
      DATA (JTAB(I),I =   2,  2)/
     :  3/
      DATA (NTAB(I),I =   4,  6)/
     :  1,0, 2/
*
*   p, d-   (j = 3/2)
*
      DATA (ITAB(I),I =   3,  4)/
     :  1,2/
      DATA (JTAB(I),I =   3,  4)/
     :    6,  9/
      DATA (NTAB(I),I =   7, 15)/
     :  1,0, 4,
     :  0,0, 1,
     :  2,0, 5/
*
*  d, f-   (j = 5/2)
*
      DATA (ITAB(I),I =   5,  7)/
     :  1,3,3/
      DATA (JTAB(I),I =   5,  7)/
     :   15, 18, 27/
      DATA (NTAB(I),I =  16, 36)/
     :  1,0, 6,
     :  0,0, 1,
     :  2,0, 5, 2,0, 9,
     :  1,0, 6,
     :  3,0, 4, 3,0,10/
*
*   f, g-   (j = 7/2)
*
      DATA (ITAB(I),I =   8, 11)/
     :  1,4,6,8/
      DATA (JTAB(I),I =   8, 11)/
     :   36, 39, 51, 69/
      DATA (NTAB(I),I =  37, 93)/
     :  1,0, 8,
     :  0,0, 1,
     :  2,0, 5, 2,0, 9, 2,0,13,
     :  1,0, 8,
     :  3,0, 4, 3,0, 6, 3,0,10, 3,0,12, 3,0,16,
     :  0,0, 1,
     :  2,0, 5, 2,0, 9, 2,0,13,
     :  4,0, 5, 4,0, 9, 4,0,11, 4,0,17/
*
*   g, h-   (j = 9/2)
*
      DATA (ITAB(I),I =  12, 16)/
     :  1,5,10,18,20/
      DATA (JTAB(I),I =  12, 16)/
     :   93, 96, 111,141,195/
      DATA (NTAB(I),I =  94,255)/
     :  1,0,10,
     :  0,0, 1,
     :  2,0, 5, 2,0, 9, 2,0,13, 2,0,17,
     :  1,0,10,
     :  3,0, 4, 3,0, 6, 3,0, 8, 3,0,10, 3,0,12, 3,0,14, 3,0,16, 3,0,18,
     :  3,0,22,
     :  0,0, 1,
     :  2,0,5, 2,0,9, 2,0,13, 2,0,17,
     :  4,0, 1, 4,0, 5, 4,0, 7, 4,0, 9, 4,1, 9, 4,0,11, 4,0,13, 4,1,13,
     :  4,0,15, 4,0,17, 4,0,19, 4,0,21, 4,0,25,
     :  1,0,10,
     :  3,0, 4, 3,0, 6, 3,0, 8, 3,0,10, 3,0,12, 3,0,14, 3,0,16, 3,0,18,
     :  3,0,22,
     :  5,0, 2, 5,0, 6, 5,0, 8, 5,0,10, 5,0,12, 5,0,14, 5,0,16, 5,0,18,
     :  5,0,20, 5,0,26/
*
*   h, i-   (j = 11/2)
*
*   First two rows only
*
      DATA (ITAB(I),I =  17, 18)/
     :  1,6/
      DATA (JTAB(I),I =  17, 19)/
     :  255,258,277/
      DATA (NTAB(I),I = 256,276)/
     :  1,0,12,
     :  0,0, 1,
     :  2,0, 5, 2,0, 9, 2,0,13, 2,0,17, 2,0,21/
*
*   i, k-   (j = 13/2)
*
*   First two rows only
*
      DATA (ITAB(I),I =  23, 24)/
     :  1,7/
      DATA (JTAB(I),I =  23, 25)/
     :  276,279,301/
      DATA (NTAB(I),I = 277,300)/
     :  1,0,14,
     :  0,0, 1,
     :  2,0, 5, 2,0, 9, 2,0,13, 2,0,17, 2,0,21, 2,0,25/
*
*   k, l-   (j = 15/2)
*
*   First two rows only
*
      DATA (ITAB(I),I =  30, 31)/
     :  1,8/
      DATA (JTAB(I),I =  30, 32)/
     :  300,303,328/
      DATA (NTAB(I),I = 301,327)/
     :  1,0,16,
     :  0,0, 1,
     :  2,0, 5, 2,0, 9, 2,0,13, 2,0,17, 2,0,21, 2,0,25, 2,0,29/
*
      END
