      MODULE convrt_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  10:46:53   2/14/04  
      SUBROUTINE convrt (INTNUM, CNUM, LENTH) 
      INTEGER, INTENT(IN) :: INTNUM 
      CHARACTER (LEN = *), INTENT(INOUT) :: CNUM 
      INTEGER, INTENT(OUT) :: LENTH 
!...This routine performs I/O.
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
      MODULE convrt_DOUBLE_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  10:46:53   2/14/04  
      SUBROUTINE convrt_double (INTNUM, CNUM, LENTH) 
      INTEGER, INTENT(IN) :: INTNUM 
      CHARACTER (LEN = *), INTENT(INOUT) :: CNUM 
      INTEGER, INTENT(OUT) :: LENTH 
!...This routine performs I/O.
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
      MODULE convrtnl_I
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90 4.3E  10:05:49  11/16/01  
      SUBROUTINE convrtnl (I, J, STRING_NL)
      INTEGER, INTENT(IN) :: I
      INTEGER, INTENT(IN) :: J
      CHARACTER (LEN = 4), INTENT(INOUT) :: STRING_NL
!VAST...Calls: 
      END SUBROUTINE
      END INTERFACE
      END MODULE
      MODULE dracah_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  07:27:24   2/14/04  
      SUBROUTINE dracah (I, J, K, L, M, N, RAC) 
      USE vast_kind_param,ONLY: DOUBLE 
      INTEGER MFACT 
      PARAMETER (MFACT = 500) 
      INTEGER, INTENT(IN) :: I 
      INTEGER, INTENT(IN) :: J 
      INTEGER, INTENT(IN) :: K 
      INTEGER, INTENT(IN) :: L 
      INTEGER, INTENT(IN) :: M 
      INTEGER, INTENT(IN) :: N 
      REAL(DOUBLE), INTENT(OUT) :: RAC 
!VAST.../FACTS/ GAM(IN)
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
      MODULE factt_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  10:47:26   2/14/04  
      SUBROUTINE factt 
!VAST.../FACTS/ GAM(INOUT)
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
      MODULE getmixblock_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  18:32:57   1/ 6/07  
      SUBROUTINE getmixblock (NAME, NCI) 
      CHARACTER (LEN = 24), INTENT(IN) :: NAME 
      INTEGER, INTENT(IN) :: NCI 
!VAST.../DEF1/ NELEC(INOUT)
!VAST.../EIGVAL/ EAV(OUT), PNEVAL(INOUT)
!VAST.../EIGVEC/ PNEVEC(INOUT)
!VAST.../ORB2/ NCF(OUT), NW(INOUT)
!VAST.../PRNT/ NVEC(OUT), PNIVEC(INOUT)
!VAST.../SYMA/ PIATJP(INOUT), PIASPA(INOUT)
!VAST.../IOUNIT/ ISTDE(IN)
!VAST...Calls: OPENFL, ALLOC, IVEC
!...This routine performs I/O.
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
      MODULE getyn_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  10:48:16   2/14/04  
      LOGICAL FUNCTION getyn ( ) 
!VAST.../IOUNIT/ ISTDI(IN), ISTDE(IN)
!...This routine performs I/O.
      END FUNCTION  
      END INTERFACE 
      END MODULE 
      MODULE ichop_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  10:48:25   2/14/04  
      INTEGER FUNCTION ichop (ISUBSH, ICSF) 
      INTEGER NNNW 
      PARAMETER(NNNW=120) 
      INTEGER, INTENT(IN) :: ISUBSH 
      INTEGER :: ICSF 
!VAST.../ORB5/ NKJ(IN)
!VAST...Calls: IQ
      END FUNCTION  
      END INTERFACE 
      END MODULE 
      MODULE idigit_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  10:48:52   2/14/04  
      INTEGER FUNCTION idigit (CST) 
      CHARACTER (LEN = 1), INTENT(IN) :: CST 
      END FUNCTION  
      END INTERFACE 
      END MODULE 
      MODULE iq_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  10:48:38   2/14/04  
      INTEGER FUNCTION iq (ISUBSH, ICSF) 
      INTEGER, INTENT(IN) :: ISUBSH 
      INTEGER :: ICSF 
!VAST...Calls: IQA
      END FUNCTION  
      END INTERFACE 
      END MODULE 
      MODULE ispar_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  10:48:41   2/14/04  
      INTEGER FUNCTION ispar (ICSF) 
      INTEGER, INTENT(IN) :: ICSF 
!VAST.../ORB2/ NCF(IN)
!VAST.../IOUNIT/ ISTDE(IN)
!VAST...Calls: JCUPA
!...This routine performs I/O.
      END FUNCTION  
      END INTERFACE 
      END MODULE 
      MODULE itjpo_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  10:48:45   2/14/04  
      INTEGER FUNCTION itjpo (ICSF) 
      INTEGER, INTENT(IN) :: ICSF 
!VAST.../ORB2/ NCF(IN)
!VAST.../IOUNIT/ ISTDE(IN)
!VAST...Calls: JCUPA
!...This routine performs I/O.
      END FUNCTION  
      END INTERFACE 
      END MODULE 
      MODULE ittk_I
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90 4.3E  09:56:17  11/16/01  
      INTEGER FUNCTION ittk (I, J, K)
      INTEGER, INTENT(IN) :: I
      INTEGER, INTENT(IN) :: J
      INTEGER, INTENT(IN) :: K
      END FUNCTION
      END INTERFACE
      END MODULE
      MODULE ixjtik_I
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90 4.3E  10:05:49  11/16/01  
      INTEGER FUNCTION ixjtik (I, J, K, L, M, N)
      INTEGER :: I
      INTEGER :: J
      INTEGER :: K
      INTEGER :: L
      INTEGER :: M
      INTEGER :: N
!VAST...Calls: ITTK
      END FUNCTION
      END INTERFACE
      END MODULE
      MODULE jcup_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  08:16:51   2/21/04  
      INTEGER FUNCTION jcup (LOC, ICSF) 
      INTEGER, INTENT(IN) :: LOC 
      INTEGER, INTENT(IN) :: ICSF 
!VAST.../ORB2/ NCF(IN), NW(IN)
!VAST.../IOUNIT/ ISTDE(IN)
!VAST...Calls: JCUPA
!...This routine performs I/O.
      END FUNCTION  
      END INTERFACE 
      END MODULE 
      MODULE jqs_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  10:48:50   2/14/04  
      INTEGER FUNCTION jqs (IWHICH, ISUBSH, ICSF) 
      INTEGER :: IWHICH 
      INTEGER, INTENT(IN) :: ISUBSH 
      INTEGER :: ICSF 
!VAST...Calls: JQSA
      END FUNCTION  
      END INTERFACE 
      END MODULE 
      MODULE lodcsl_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  13:07:22   2/14/04  
      SUBROUTINE lodcsl (NCORE) 
      INTEGER NNNW 
      PARAMETER(NNNW=120) 
      INTEGER NW2 
      PARAMETER(NW2=240) 
      INTEGER, INTENT(OUT) :: NCORE 
!VAST.../DEBUGA/ LDBPA(IN)
!VAST.../DEF1/ NELEC(OUT)
!VAST.../ORB2/ NCF(OUT), NW(OUT), PNTRIQ(INOUT)
!VAST.../ORB4/ NP(INOUT), NAK(INOUT)
!VAST.../ORB5/ NKL(INOUT), NKJ(INOUT)
!VAST.../ORB10/ NH(INOUT)
!VAST.../STAT/ PNTJQS(INOUT), PNJCUP(INOUT)
!VAST.../TERMS/ JTAB(IN), NTAB(IN)
!VAST.../IOUNIT/ ISTDE(IN)
!VAST.../BLK/ NBLOCK(OUT), NCFBLK(OUT)
!VAST...Calls: PRSRSL, CONVRT, ALLOC, PRSRCN, PARSJL, RALC2D
!VAST...Calls: PACK, IQ, JQS, JCUP, ITJPO, ISPAR
!...This routine performs I/O.
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
      MODULE lodiso_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  10:49:53   2/14/04  
      SUBROUTINE lodiso 
!VAST.../DEF1/ EMN(OUT), Z(INOUT)
!VAST.../DEF11/ FMTOAU(IN), AUMAMU(IN)
!VAST.../NPAR/ PARM(OUT), NPARM(OUT)
!VAST.../NSMDAT/ SQN(INOUT), DMOMNM(INOUT), QMOMB(INOUT)
!...This routine performs I/O.
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
      MODULE lval_I
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90 4.3E  21:49:59  11/14/01
      INTEGER FUNCTION lval (SYMBOL)
      CHARACTER (LEN = 1), INTENT(IN) :: SYMBOL
      END FUNCTION
      END INTERFACE
      END MODULE
      MODULE nine0_I
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90 4.3E  10:07:58  11/16/01  
      SUBROUTINE nine0 (J1, J2, J3, L1, L2, L3, K1, K2, K3, AA)
      USE vast_kind_param,ONLY: DOUBLE
      INTEGER, INTENT(IN) :: J1
      INTEGER, INTENT(IN) :: J2
      INTEGER, INTENT(IN) :: J3
      INTEGER, INTENT(IN) :: L1
      INTEGER, INTENT(IN) :: L2
      INTEGER, INTENT(IN) :: L3
      INTEGER, INTENT(IN) :: K1
      INTEGER, INTENT(IN) :: K2
      INTEGER, INTENT(IN) :: K3
      REAL(DOUBLE), INTENT(OUT) :: AA
!VAST.../CONSTS/ ZERO(IN), ONE(IN)
!VAST...Calls: SIXJ
      END SUBROUTINE
      END INTERFACE
      END MODULE
      MODULE nine_I
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90 4.3E  10:07:58  11/16/01  
      SUBROUTINE nine (J1, J2, J3, L1, L2, L3, K1, K2, K3, I, IN, AA)
      USE vast_kind_param,ONLY: DOUBLE
      INTEGER, INTENT(IN) :: J1
      INTEGER, INTENT(IN) :: J2
      INTEGER, INTENT(IN) :: J3
      INTEGER, INTENT(IN) :: L1
      INTEGER, INTENT(IN) :: L2
      INTEGER, INTENT(IN) :: L3
      INTEGER, INTENT(IN) :: K1
      INTEGER, INTENT(IN) :: K2
      INTEGER, INTENT(IN) :: K3
      INTEGER, INTENT(IN) :: I
      INTEGER, INTENT(OUT) :: IN
      REAL(DOUBLE), INTENT(INOUT) :: AA
!VAST.../CONSTS/ ZERO(IN)
!VAST...Calls: ITTK, NINE0, SIXJ
      END SUBROUTINE
      END INTERFACE
      END MODULE
      MODULE openfl_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  10:50:04   2/14/04  
      SUBROUTINE openfl (NFILE, FILNAM, RFORM, RSTAT, IERR) 
      INTEGER, INTENT(IN) :: NFILE 
      CHARACTER (LEN = *), INTENT(IN) :: FILNAM 
      CHARACTER (LEN = *), INTENT(IN) :: RFORM 
      CHARACTER (LEN = *), INTENT(IN) :: RSTAT 
      INTEGER, INTENT(OUT) :: IERR 
!VAST.../IOUNIT/ ISTDE(IN)
!...This routine performs I/O.
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
      MODULE packLS_I
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90 4.3E  21:49:59  11/14/01
      SUBROUTINE packLS (M, EL, Q, COUPLE, STR)
      INTEGER, INTENT(IN) :: M
      CHARACTER (LEN = 3), DIMENSION(*), INTENT(IN) :: EL
      INTEGER, DIMENSION(*), INTENT(IN) :: Q
      CHARACTER (LEN = 3), DIMENSION(*), INTENT(IN) :: COUPLE
      CHARACTER (LEN = 64), INTENT(OUT) :: STR
!VAST...Calls: LVAL
      END SUBROUTINE
      END INTERFACE
      END MODULE
      MODULE pack_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  10:50:08   2/14/04  
      SUBROUTINE pack (IUNPKD, ISUBSH, IPACKD) 
      USE vast_kind_param,ONLY: BYTE 
      INTEGER, INTENT(IN) :: IUNPKD 
      INTEGER, INTENT(IN) :: ISUBSH 
      INTEGER(BYTE), DIMENSION(*), INTENT(INOUT) :: IPACKD 
!VAST.../IOUNIT/ ISTDE(IN)
!...This routine performs I/O.
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
      MODULE parsjl_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  10:50:09   2/14/04  
      SUBROUTINE parsjl (MODE, NCORE, RECORD, LOC, JX, NJX, IERR) 
      INTEGER, INTENT(IN) :: MODE 
      INTEGER, INTENT(IN) :: NCORE 
      CHARACTER (LEN = 256), INTENT(IN) :: RECORD 
      INTEGER, INTENT(IN) :: LOC 
      INTEGER, DIMENSION(*), INTENT(OUT) :: JX 
      INTEGER, INTENT(OUT) :: NJX 
      INTEGER, INTENT(OUT) :: IERR 
!VAST.../ORB2/ NW(IN)
!VAST...Calls: CONVRT
!...This routine performs I/O.
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
      MODULE prsrcn_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  10:50:15   2/14/04  
      SUBROUTINE prsrcn (RECORD, NCORE, IOCCS, IERR) 
      INTEGER NNNW 
      PARAMETER(NNNW=120) 
      CHARACTER (LEN = 256), INTENT(IN) :: RECORD 
      INTEGER, INTENT(IN) :: NCORE 
      INTEGER, DIMENSION(NNNW), INTENT(OUT) :: IOCCS 
      INTEGER, INTENT(OUT) :: IERR 
!VAST.../ORB2/ NW(IN)
!VAST.../ORB4/ NP(IN)
!VAST.../ORB5/ NKJ(IN)
!VAST.../ORB10/ NH(IN)
!...This routine performs I/O.
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
      MODULE prsrsl_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  10:50:16   2/14/04  
      SUBROUTINE prsrsl (NFILE, ID) 
      INTEGER NNNW 
      PARAMETER(NNNW=120) 
      INTEGER, INTENT(IN) :: NFILE 
      INTEGER, INTENT(IN) :: ID 
!VAST.../IOUNIT/ ISTDE(IN)
!VAST.../ORB2/ NW(INOUT)
!VAST.../ORB4/ NP(INOUT), NAK(OUT)
!VAST.../ORB5/ NKL(OUT), NKJ(OUT)
!VAST.../ORB10/ NH(OUT)
!VAST...Calls: CONVRT
!...This routine performs I/O.
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
      MODULE setcsla_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  10:50:33   2/14/04  
      SUBROUTINE setcsla (NAME, NCORE) 
      CHARACTER (LEN = 24), INTENT(IN) :: NAME 
      INTEGER :: NCORE 
!VAST.../IOUNIT/ ISTDE(IN)
!VAST...Calls: OPENFL, LODCSL
!...This routine performs I/O.
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
      MODULE setiso_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  10:50:35   2/14/04  
      SUBROUTINE setiso (FNAME) 
      CHARACTER (LEN = *), INTENT(INOUT) :: FNAME 
!VAST.../IOUNIT/ ISTDE(IN)
!VAST...Calls: OPENFL, LODISO
!...This routine performs I/O.
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
      MODULE sixj1_I
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90 4.3E  10:06:26  11/16/01  
      SUBROUTINE sixj1 (I, J, K, L, M, ITIK, SI)
      USE vast_kind_param,ONLY: DOUBLE
      INTEGER, INTENT(IN) :: I
      INTEGER, INTENT(IN) :: J
      INTEGER, INTENT(IN) :: K
      INTEGER, INTENT(IN) :: L
      INTEGER, INTENT(IN) :: M
      INTEGER, INTENT(IN) :: ITIK
      REAL(DOUBLE), INTENT(OUT) :: SI
!VAST.../CONSTS/ ZERO(IN), HALF(IN), ONE(IN), TWO(IN), THREE(IN)
!VAST...Calls: IXJTIK
      END SUBROUTINE
      END INTERFACE
      END MODULE
      MODULE sixj2_I
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90 4.3E  10:06:26  11/16/01  
      SUBROUTINE sixj2 (J, K, L, M, N, ITIK, SI)
      USE vast_kind_param,ONLY: DOUBLE
      INTEGER, INTENT(IN) :: J
      INTEGER, INTENT(IN) :: K
      INTEGER, INTENT(IN) :: L
      INTEGER, INTENT(IN) :: M
      INTEGER, INTENT(IN) :: N
      INTEGER, INTENT(IN) :: ITIK
      REAL(DOUBLE), INTENT(OUT) :: SI
!VAST.../CONSTS/ ZERO(IN), HALF(IN), ONE(IN), TWO(IN), THREE(IN) 
!VAST.../CONSTS/ FOUR(IN)
!VAST...Calls: IXJTIK
      END SUBROUTINE
      END INTERFACE
      END MODULE
      MODULE sixj35_I
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90 4.3E  10:06:26  11/16/01  
      SUBROUTINE sixj35 (J, K, L, M, N, ITIK, SI)
      USE vast_kind_param,ONLY: DOUBLE
      INTEGER, INTENT(IN) :: J
      INTEGER, INTENT(IN) :: K
      INTEGER, INTENT(IN) :: L
      INTEGER, INTENT(IN) :: M
      INTEGER, INTENT(IN) :: N
      INTEGER, INTENT(IN) :: ITIK
      REAL(DOUBLE), INTENT(OUT) :: SI
!VAST.../CONSTS/ ZERO(IN), ONE(IN), TWO(IN), THREE(IN), FOUR(IN)
!VAST...Calls: IXJTIK
      END SUBROUTINE
      END INTERFACE
      END MODULE
      MODULE sixj3_I
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90 4.3E  10:06:26  11/16/01  
      SUBROUTINE sixj3 (J, K, L, M, N, ITIK, SI)
      USE vast_kind_param,ONLY: DOUBLE
      INTEGER, INTENT(IN) :: J
      INTEGER, INTENT(IN) :: K
      INTEGER, INTENT(IN) :: L
      INTEGER, INTENT(IN) :: M
      INTEGER, INTENT(IN) :: N
      INTEGER, INTENT(IN) :: ITIK
      REAL(DOUBLE), INTENT(OUT) :: SI
!VAST.../CONSTS/ ZERO(IN), ONE(IN), TWO(IN), THREE(IN), FOUR(IN) 
!VAST.../CONSTS/ SEVEN(IN)
!VAST...Calls: IXJTIK
      END SUBROUTINE
      END INTERFACE
      END MODULE
      MODULE sixj4_I
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90 4.3E  10:06:26  11/16/01  
      SUBROUTINE sixj4 (JC, JE, JD, JB, JF, ITIK, SI)
      USE vast_kind_param,ONLY: DOUBLE
      INTEGER, INTENT(IN) :: JC
      INTEGER, INTENT(IN) :: JE
      INTEGER, INTENT(IN) :: JD
      INTEGER, INTENT(IN) :: JB
      INTEGER, INTENT(IN) :: JF
      INTEGER, INTENT(IN) :: ITIK
      REAL(DOUBLE), INTENT(OUT) :: SI
!VAST.../CONSTS/ ZERO(IN), HALF(IN), ONE(IN), TWO(IN), THREE(IN)
!VAST.../CONSTS/ EPS(IN)
!VAST...Calls: IXJTIK, GRACAH1, SIXJ2, SIXJ3
      END SUBROUTINE
      END INTERFACE
      END MODULE
      MODULE sixj5_I
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90 4.3E  10:06:26  11/16/01  
      SUBROUTINE sixj5 (J, K, L, M, N, ITIK, SI)
      USE vast_kind_param,ONLY: DOUBLE
      INTEGER, INTENT(IN) :: J
      INTEGER, INTENT(IN) :: K
      INTEGER, INTENT(IN) :: L
      INTEGER, INTENT(IN) :: M
      INTEGER, INTENT(IN) :: N
      INTEGER, INTENT(IN) :: ITIK
      REAL(DOUBLE), INTENT(OUT) :: SI
!VAST.../CONSTS/ ZERO(IN), ONE(IN), TWO(IN)
!VAST...Calls: IXJTIK
      END SUBROUTINE
      END INTERFACE
      END MODULE
      MODULE sixj_I
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90 4.3E  10:06:26  11/16/01
      SUBROUTINE sixj (I, J, K, L, M, N, ITIK, SI)
      USE vast_kind_param,ONLY: DOUBLE
      INTEGER, INTENT(IN) :: I
      INTEGER, INTENT(IN) :: J
      INTEGER, INTENT(IN) :: K
      INTEGER, INTENT(IN) :: L
      INTEGER, INTENT(IN) :: M
      INTEGER, INTENT(IN) :: N
      INTEGER, INTENT(IN) :: ITIK
      REAL(DOUBLE), INTENT(OUT) :: SI
!VAST.../CONSTS/ ZERO(IN), ONE(IN)
!VAST...Calls: IXJTIK, SIXJ5, SIXJ1, SIXJ35, SIXJ2, GRACAH1
!VAST...Calls: SIXJ3, SIXJ4
      END SUBROUTINE
      END INTERFACE
      END MODULE