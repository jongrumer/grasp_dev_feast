************************************************************************
*                                                                      *
      SUBROUTINE NJGRAF (RECUP,IGEN,FAIL)
*                                                                      *
*   Program to calculate a general recoupling coefficient. This ver-   *
*   sion is slightly modified  (Colin T Johnson, Oxford University).   *
*   The changes are as follows:                                        *
*                                                                      *
*      1. PARAMETER IGEN has been included in the argument list:       *
*            IGEN =  0  normal call to NJGRAF                          *
*            IGEN = -1  GENSUM is not called                           *
*      2. The contents of COMMON blocks /ARGU/ and /SUMARG/ are used   *
*         by  GENSUM  to calculate the recoupling coefficient. These   *
*         COMMON  blocks  have  been removed from GENSUM. Their con-   *
*         tents are passed to  GENSUM  through the argument list in-   *
*         stead, so that NJGRAF can be called to set up formulae for   *
*         both the direct and exchange cases in COR and BREIT.         *
*      3. Extra  dimension  tests  have  been  included  in routines   *
*         NJGRAF, PRINTJ, SPRATE, VAR and ZERO. These are  discussed   *
*         below.                                                       *
*      4. An  extra routine  RDIAG  has been introduced to remove an   *
*         extended  DO loop from GENSUM, to conform with the FORTRAN   *
*         77 standard.                                                 *
*                                                                      *
*                                                                      *
*   Description of some COMMON blocks. A full discussion is given in   *
*   the NJGRAF program description (Bar-Shalom and Klapisch op cit).   *
*                                                                      *
*      COMMON block COUPLE                                             *
*                                                                      *
*         M                The total number of angular momentum val-   *
*                          ues in the initial and final states         *
*         N                The number of basic angular momentum val-   *
*                          ues that are coupled                        *
*         J1(I),           The angular momentum values stored as 2J+1  *
*            I = 1,M                                                   *
*         J2(I,J),         The position in the J1 array of the init-   *
*            I = 1,(N-1),  ial state triads                            *
*            J = 1,3                                                   *
*         J3(I,J),         The position in the J1 array of the final   *
*            I = 1,(N-1),  state triads                                *
*            J = 1,3                                                   *
*         FREE(I),         If FREE(I) = .TRUE., no reference is made   *
*            I = 1,M       to the value of J1(I) when establishing a   *
*                          formula in  NJGRAF .  GENSUM  may then be   *
*                          called  for  repeated  occurences of this   *
*                          formula  with  differing values of J1(I).   *
*                          If J1(I) does  not  vary between calls to   *
*                          GENSUM then FREE(I) should be set .FALSE.   *
*                          so that zero branches  can be removed be-   *
*                          fore the formula is established.            *
*                                                                      *
*      COMMON block DEBUG                                              *
*                                                                      *
*         IBUG1            Not used                                    *
*         IBUG2            Not used                                    *
*         IBUG3            Debug prints in NJGRAF and GENSUM if 1      *
*         IBUG4            Not used                                    *
*         IBUG5            Not used                                    *
*         IBUG6            Not used                                    *
*                                                                      *
*      COMMON block ARGU                                               *
*                                                                      *
*         J6C              The number of elements in the K6 array      *
*         J7C              The number of elements in the K7 array      *
*         J8C              The number of elements in the K8 array      *
*         J9C              The number of elements in the K9 array      *
*         JWC              The number of columns in the KW array       *
*         J6(I),           Each entry corresponds to a  factor  SQRT   *
*            I = 1,J6C     (2J+1) in RECUP. The value  of  J6  GIVES   *
*                          position in  J1  array where  J  value is   *
*                          found                                       *
*         J7(I),           Each entry corresponds to factor  (-1)**J   *
*            I = 1,J7C     in RECUP                                    *
*         J8(I),           Each entry corresponds to a factor (-1)**   *
*            I = 1,J8C     (2J) in RECUP                               *
*         J9(I),           Each entry corresponds to a factor (2J+1)   *
*            I = 1,J9C     **(-0.5) in RECUP                           *
*         KW(I,J),         Each column corresponds to a Racah coeff-   *
*            I = 1,6,      icient in RECUP                             *
*            J = 1,JWC                                                 *
*         JDEL             The number of delta functions               *
*         LDEL(I,J),       The arguments of the delta functions        *
*              J = 1,2                                                 *
*         SUMVAR(I)        .TRUE. for ang. mom. I (a summation vari-   *
*                          able                                        *
*         MP               The index of the last variable              *
*                                                                      *
*   The arrays  J6, J7, J8, J9 and  KW, Are evaluated by NJGRAF. The   *
*   summation over the variables in  J6, J7, J8, J9 and  KW, and the   *
*   evaluation of RECUP is carried out in GENSUM. GENSUM  can be re-   *
*   entered directly to evaluate different  recoupling  coefficients   *
*   with the same structure  by just  altering the numbers in the J1   *
*   array.                                                             *
*                                                                      *
*   This is the main program. It handles all the analysis of the re-   *
*   coupling  coefficient without referring explicitly to the values   *
*   of angular  momenta  which  are in J1(J),except for zero in case   *
*   FREE = .FALSE. . Like NJSYM it  prepares arrays of arguments for   *
*   phase factors, (2*J+1) factors and  6j-coefficients to be compu-   *
*   ted in GENSUM, which can also be called separately when only the   *
*   numerical values of angular momenta change. These variable angu-   *
*   lar momenta should be declared  FREE(J)  = .TRUE. , so  that the   *
*   formula prepared for GENSUM should be correct when J1 is not ze-   *
*   ro. FAIL will be TRUE when the  recoupling  coefficient  is zero   *
*   because of unsatisfied delta or other similar causes.              *
*                                                                      *
*   This version holds the array dimensions in parameter statements.   *
*   The dimensions are labelled:                                       *
*                                                                      *
*      MANGM  : Dimension of the J1 and FREE arrays in /COUPLE/, and   *
*               the  first  dimension of the LINE and LCOL arrays in   *
*               /TREE/. Also  the  dimension  of the SUMVAR array in   *
*               /ARGU/, AND OF THE INVER array in routine SPRATE. It   *
*               is tested for  M  on entry to  NJGRAF, and for MP in   *
*               routine SPRATE.                                        *
*      MTRIAD : Dimension of the  J2 and  J3 arrays in /COUPLE/. The   *
*               dimensions of these  arrays  are checked on entry to   *
*               NJGRAF in addition  MTRIAD sets the dimension of the   *
*               JSUM6 array and the first dimension of the JSUM4 and   *
*               JSUM5  arrays in /SUMARG/. Also gives the dimensions   *
*               of some  temporary working arrays in SPRATE and GEN-   *
*               SUM. In these  cases  mtriad sets the maximum number   *
*               of summation variables  in any particular sum, which   *
*               is tested in SPRATE.                                   *
*      M2TRD  : (=2*MTRIAD) Dimension of the J23 ,  ARROW  and  TABS   *
*               arrays in /TREE/. Also  the  dimension of the npoint   *
*               array in /GRAPH/.                                      *
*      M4TRD  : (=4*MTRIAD) Dimension of the  JDIAG,  ARR, IL and IH   *
*               arrays in /GRAPH/, and of the IAL array in /BUILD/.    *
*      M3MNGM : Dimension of the J6 array in /ARGU/, tested in SPRATE  *
*               Dimension of the J7 array in /ARGU/, tested in SPRATE  *
*               Dimension of the J8 array in /ARGU/, tested in SPRATE  *
*      MANGMP : Dimension of the J9 array in /ARGU/, tested in SPRATE  *
*               MANGMP also sets the dimension of the J6P,  J7P, J8P   *
*               and J9P arrays in /SUMARG/, And of the JNS  array in   *
*               routine VAR. The dimension of the JNS array is  tes-   *
*               ted in VAR.                                            *
*      M6J    : Dimension of the JW(or KW) and LDEL arrays in /ARGU/,  *
*               and of the JWORD and INV6J arrays in /SUMARG/.  Also   *
*               the second dimension of the  JSUM4 and  JSUM5 arrays   *
*               in /SUMARG/. In addition it  gives the dimensions of   *
*               a  number of  temporary  working  arrays in routines   *
*               SPRATE and GENSUM. M6J is tested in SPRATE.            *
*      MFACT  : The dimension of the factorial array GAM in /FACTS /.  *
*      MSUM   : Dimension of the NBJ, NB6J, K6CP, K7CP, K8CP and K9CP  *
*               arrays in /SUMARG/. MSUM is the  maximum  number  of   *
*               sums allowed, and is tested in routine SPRATE.         *
*      MTAB   : The dimension of the JTAB array in  routine  PRINTJ.   *
*               MTAB is tested in PRINTJ.                              *
*      MZERO  : Dimension of the JZERO array in /ZER/. MZERO is tes-   *
*               ted in routine ZERO.                                   *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      CHARACTER*6 NAME,NAMSUB
      LOGICAL FAIL,FIND,TABS,CUT,FREE,SUMVAR
      INTEGER ARROW,ARR,TAB1
*
      PARAMETER (MTRIAD = 12, M2TRD = 2*MTRIAD, M4TRD = 4*MTRIAD,
     :           MANGM = 60,  M3MNGM = 3*MANGM, MANGMP = 2*(MANGM/3),
     :           MFACT = 500, M6J = 20, MSUM = 10)
*
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :            J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),
     :            MP
     :      /CONS/ZRO,HALF,TENTH,ONE,TWO,THREE,TEN
     :      /COUPLE/M,N,J1(MANGM),J2(MTRIAD,3),J3(MTRIAD,3),FREE(MANGM)
     :      /CUTDIG/CUT
     :      /DEBUG/IBUG1,IBUG2,IBUG3,IBUG4,IBUG5,IBUG6
     :      /FACTS /GAM(MFACT)
     :      /GRAPH/JDIAG(M4TRD,3),ARR(M4TRD,3),TAB1(MANGM,2),IL(M4TRD),
     :             IH(M4TRD),NPOINT(M2TRD),NBNODE,IFIRST,ILAST,IPARTS,
     :             IPARTL,NPART,ICROSS,NFREE,ITFREE(M6J),NFIN,NC
     :      /NAM/NAMSUB
     :      /SUMARG/J6P(MANGMP),J7P(MANGMP),J8P(MANGMP),J9P(MANGMP),
     :              JWORD(6,M6J),NLSUM,NBJ(MSUM),NB6J(MSUM),K6CP(MSUM),
     :              K7CP(MSUM),K8CP(MSUM),K9CP(MSUM),JSUM6(MTRIAD),
     :              JSUM4(MTRIAD,M6J),JSUM5(MTRIAD,M6J),INV6J(M6J)
     :      /TREE/J23(M2TRD,3),ARROW(M2TRD,3),LINE(MANGM,2),
     :            LCOL(MANGM,2),TABS(M2TRD),NBTR
*
      DATA NAME/'NJGRAF'/
*
*   Debug printout
*
      IF (IBUG3 .EQ. 1) THEN
         WRITE (99,300)
         WRITE (99,301) M,N-1
         WRITE (99,302)
         WRITE (99,303) (J1(I),I = 1,M)
         WRITE (99,304) (FREE(I),I = 1,M)
         WRITE (99,305)
         WRITE (99,306) ((J2(I,J),J = 1,3),
     :                    (J3(I,J),J = 1,3),I = 1,N-1)
      ENDIF
*
*   Test the dimension of the J1 array
*
      IF (M+1 .GT. MANGM) THEN
         WRITE (*,307)
         WRITE (*,308) M+1,MANGM
         STOP
      ENDIF
*
*   Test the dimensions of the J2 and J3 arrays
*
      IF (N-1 .GT. MTRIAD) THEN
         WRITE (*,307)
         WRITE (*,309) N-1,MTRIAD
         STOP
      ENDIF
*
*   Initializations
*
      DO 502 I = N,MTRIAD
         DO 501 J = 1,3
            J2(I,J) = 0
            J3(I,J) = 0
  501    CONTINUE
  502 CONTINUE
*
      FAIL = .FALSE.
      J6C = 0
      J7C = 0
      J8C = 0
      J9C = 0
      JWC = 0
      JDEL = 0
      CALL SETDM
      NFIN = 0
      CUT = .FALSE.
*
*   Building up of the unstructured graph
*
      CALL SETTAB (FAIL)
*
*   Exit with RECUP set to zero if any angular momentum is
*   impossible
*
      M = M+1
      IF (FAIL) GOTO 7
*
      M = M-1
*
*   Locate and eliminate any zero angular momentum; simplify the
*   graph
*
      JF = 0
      JF1 = 0
      CALL ZERO (JF1,JF,FAIL)
      IF (FAIL) GOTO 7
*
      MP = M
      IF (NBTR .EQ. 0) GOTO 6
      JUMP = 1
*
    1 CALL CHKLP1 (FAIL)
      IF (FAIL) GOTO 7
*
*   Build a flat diagram out of the unstructured graph; several flat
*   diagrams may constitute the original graph, in which case there
*   are possible cuts; the flat diagrams will have free ends if cut
*
      CALL DIAGRM (JUMP)
      NFIN = MAX (0,NFREE-2)
*
      IF (NFIN .NE. 0) THEN
         JUMP = 3
*
*   Handling of free ends if a cut was found
*
         CALL CUTNL (FAIL)
         IF (FAIL) GOTO 7
      ELSE
         JUMP = 2
         IF (NFREE .EQ. 1) THEN
            CALL CUT1L (FAIL)
            IF (FAIL) GOTO 7
         ELSEIF (NFREE .GT. 1) THEN
            CALL CUT2L (FAIL)
            IF (FAIL) GOTO 7
         ENDIF
      ENDIF
*
      NBTR = NBTR+NFIN
      IF (NBTR .NE. 0) CUT = .TRUE.
*
*   Analysis of the flat diagram.
*   Closed circuits of increasing order NC are searched, analysed,
*   and taken out of the flat diagram, thus reducing the number of
*   nodes, NBNODE.
*
      NC = 0
   10 NC = NC+1
      CALL SEARCH (FIND)
      IF (.NOT. FIND) GOTO 10
      NCP = NC-2
      JPOL = 0
      IF ((M .EQ. MP) .AND. (NC.GT.3)) CALL SETDM
      IF (IPARTL .GT. 2) CALL POLYGN (JPOL)
      GOTO (11,12,13,14),NC
   11 CALL LOLPOP (FAIL)
      IF (FAIL) GOTO 7
      GOTO 15
   12 CALL BUBBLE (JPOL,FAIL)
      IF (FAIL) GOTO 7
      GOTO 15
   13 CALL TRIANG (FAIL)
      IF (FAIL) GOTO 7
      GOTO 15
   14 CALL SQUARE
   15 NBNODE = NBNODE-2
      IF (NBNODE .EQ. 0) GOTO 9
      IFIRST = IH(1)
      ILAST = IH(NBNODE)
*
*   PRINTJ is an all purpose printing SUBROUTINE called from many
*   places
*
      CALL PRINTJ (NAMSUB,8)
      IF (NBNODE .EQ. NFIN) GOTO 9
      NC = NCP
*
*   Proceed to other circuits of order NC-1
*
      GOTO 10
    9 IF (NBTR .EQ. 0) GOTO 6
      IF (JUMP .EQ. 3) CALL ORDTRI
*
*   At this stage, the flat diagram has been reduced to nodes
*   involving free ends. Proceed to build other flat diagrams
*   if necessary.
*
      GOTO 1
*
*   All parts of the original graph have been reduced.
*
    7 RECUP = ZRO
      M = M-1
      RETURN
    6 CALL PRINTJ (NAME,0)
*
*   Preparation of the results, and separation in several sums
*   if cuts have been detected, also in the flat diagram itself
*
      CALL SPRATE (M)
      M = M-1
*
*   GENSUM computes the numerical value of the recoupling
*   coefficient
*
      IF (IGEN .NE. -1) CALL GENSUM (J6C,J7C,J8C,J9C,JWC,J6,J7,J8,J9,KW,
     :                               JDEL,LDEL,SUMVAR,MP,J6P,J7P,J8P,
     :                               J9P,JWORD,NLSUM,NBJ,NB6J,K6CP,K7CP,
     :                               K8CP,K9CP,JSUM4,JSUM5,JSUM6,INV6J,
     :                               RECUP)
*
      RETURN
*
  300 FORMAT (//' ++++++++++ NJGRAF ++++++++++'/)
  301 FORMAT (' Total number of angular momenta (M) = ',1I3
     :      //' Number of triads in each of the ''left-hand'' and',
     :        ' ''right-hand'' states (N-1) = ',1I3)
  302 FORMAT (/' (2J+1)-value for each angular momentum:')
  303 FORMAT (1X,42I3)
  304 FORMAT (1X,42L3)
  305 FORMAT (/' ''Left-hand'' triads',10X,'''Right-hand'' triads')
  306 FORMAT (1X,3I3,19X,3I3)
  307 FORMAT (/' ***** Error in NJGRAF *****'/)
  308 FORMAT (' M+1 = ',1I3,', exceeds PARAMETER MANGM = ',1I3)
  309 FORMAT (' N-1 = ',1I3,', exceeds PARAMETER MTRIAD = ',1I3)
*
      END
