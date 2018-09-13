************************************************************************
*                                                                      *
      SUBROUTINE BRINT3 (IA,IB,IC,ID,NU,TEGRAL)
*                                                                      *
*   Returns integrals for the transverse photon interaction.           *
*      Integrals are stored in ordered lists. If the integral cannot   *
*   be read from a list, it is computed by calling BRINTF.             *
*                                                                      *
*   Observ that it is not possible to use more than 100 orbitals when  *
*   the Breit interaction is included. If 100 orbitals are used then   *
*   NU <= 19.                                                          *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, RALLOC.                                *
*               [RCI92]: BRINTF                                        *
*                                                                      *
*   Written by Farid A Parpia               Last update: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      LOGICAL FIRST,FOUND
*
Cww   INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)

      POINTER (PINDT1,INDTP1(*))
      POINTER (PVALT1,VALTP1(*))
      POINTER (PINDT2,INDTP2(*))
      POINTER (PVALT2,VALTP2(*))
      POINTER (PINDT3,INDTP3(*))
      POINTER (PVALT3,VALTP3(*))
      POINTER (PINDT4,INDTP4(*))
      POINTER (PVALT4,VALTP4(*))
      POINTER (PINDT5,INDTP5(*))
      POINTER (PVALT5,VALTP5(*))
      POINTER (PINDT6,INDTP6(*))
      POINTER (PVALT6,VALTP6(*))
*
      COMMON/BILST/PINDT1,PINDT2,PINDT3,PINDT4,PINDT5,PINDT6,
     :             PVALT1,PVALT2,PVALT3,PVALT4,PVALT5,PVALT6,
     :             NDTPA(6),NTPI(6),FIRST(6)
     :      /ORB2/NCF,NW,PNTRIQ
*
      KEY = NW + 1
*
*   Compute the integral label
*
      INDEX = (((NU*KEY+IA)*KEY+IB)*KEY+IC)*KEY+ID
*
      IF (.NOT. FIRST(3)) THEN
*
*   This branch is executed on all entries except the first
*
         IF (INDEX .GT. INDTP3(NTPI(3))) THEN
*
*   The index is greater than the largest stored
*
            FOUND = .FALSE.
            LOC = NTPI(3)
*
         ELSEIF (INDEX .LT. INDTP3(1)) THEN
*
*   The index is less than the smallest stored
*
            FOUND = .FALSE.
            LOC = 0
*
         ELSE
*
*   The index is within the range of the indices stored; search
*   for it in the list of indices
*
            JU = NTPI(3)
            JL = 1
    1       IF (JU-JL .GT. 1) THEN
               JM = (JU+JL)/2
               IF (INDTP3(JM) .GT. INDEX) THEN
                  JU = JM
               ELSE
                  JL = JM
               ENDIF
               GOTO 1
*
            ELSE
*
*   The range is bracketed to the extent possible
*
               IF (INDEX .EQ. INDTP3(JU)) THEN
*
                  FOUND = .TRUE.
                  LOC = JU
*
               ELSEIF (INDEX .EQ. INDTP3(JL)) THEN
*
                  FOUND = .TRUE.
                  LOC = JL
*
               ELSE
*
                  FOUND = .FALSE.
                  LOC = JL
*
               ENDIF
*
            ENDIF
*
         ENDIF
*
         IF (FOUND) THEN
*
*   Found the index in the list; return the value of the integral
*   from storage
*
            TEGRAL = VALTP3(LOC)
*
         ELSE
*
*   Index not found; compute the integral
*
            TEGRAL = BRINTF (3,IA,IB,IC,ID,NU)
*
*   Increment the integral counter
*
            NTPI(3) = NTPI(3)+1
*
*   Increase array length by half the present length if the latter
*   is inadequate to store another pair
*
*
            IF (NTPI(3) .GT. NDTPA(3)) THEN
               NEWSIZ = NDTPA(3)+NDTPA(3)/2
               CALL RALLOC (PINDT3,NDTPA(3),NEWSIZ,4)
               CALL RALLOC (PVALT3,NDTPA(3),NEWSIZ,8)
               NDTPA(3)= NEWSIZ
            ENDIF
            DO 4 I = NTPI(3),LOC+2,-1
               INDTP3(I) = INDTP3(I-1)
               VALTP3(I) = VALTP3(I-1)
    4       CONTINUE
*
*   Put the new index and value into storage
*
            INDTP3(LOC+1) = INDEX
            VALTP3(LOC+1) = TEGRAL
*
         ENDIF
*
      ELSE
*
*   This branch is executed only once per type of integral
*
         FIRST(3) = .FALSE.
*
*   Designate the initial storage for arrays INDTPx and VALTPx;
*   Array NDTPA stores the array dimensions
*
         NDTPA(3) = 10000
*
*   Compute the integral's value
*
         TEGRAL = BRINTF (3,IA,IB,IC,ID,NU)
*
*   Initialise the integral counter
*
         NTPI(3) = 1
*
*   Store the integral and its value
*
         CALL ALLOC (PINDT3,NDTPA(3),4)
         CALL ALLOC (PVALT3,NDTPA(3),8)
         INDTP3(1) = INDEX
         VALTP3(1) = TEGRAL
*
      ENDIF
*
      RETURN
      END
