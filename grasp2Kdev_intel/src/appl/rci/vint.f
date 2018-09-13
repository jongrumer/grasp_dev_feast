************************************************************************
*                                                                      *
      SUBROUTINE VINT (IA,IB,TEGRAL)
*                                                                      *
*   This routine returns Vinti integrals.                              *
*      A list of previously computed integrals is maintained. If the   *
*   required  integral is not already  available it is computed  and   *
*   stored in the list.                                                *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, RALLOC, VINTI.                         *
*                                                                      *
*   Written by Farid A Parpia             Last revision: 06 Jun 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      LOGICAL FIRST,FOUND
*
      POINTER (PVINIL,INDTEI(*))
      POINTER (PVINVL,VALTEI(*))
*
      COMMON/VINLST/NDVIN,NVINTI,PVINIL,PVINVL,FIRST
*
      include 'parameters.def'
CGG      PARAMETER (KEYORB = 121)
*
*   Compute the composite (packed) index
*
      INDEX = IA*KEYORB+IB
*
      IF (.NOT. FIRST) THEN
*
*   This branch is executed on all entries except the first
*
         IF (INDEX .GT. INDTEI(NVINTI)) THEN
*
*   The index is greater than the largest stored
*
            FOUND = .FALSE.
            LOC = NVINTI
*
         ELSEIF (INDEX .LT. INDTEI(1)) THEN
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
            JU = NVINTI
            JL = 1
*
    1       IF (JU-JL .GT. 1) THEN
*
               JM = (JU+JL)/2
               IF (INDTEI(JM) .GT. INDEX) THEN
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
               IF (INDEX .EQ. INDTEI(JU)) THEN
*
                  FOUND = .TRUE.
                  LOC = JU
*
               ELSEIF (INDEX .EQ. INDTEI(JL)) THEN
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
            TEGRAL = VALTEI(LOC)
*
         ELSE
*
*   Index not found; compute the integral
*
            TEGRAL = VINTI (IA,IB)
*
*   Increment the integral counter
*
            NVINTI = NVINTI+1
*
*   Increase array length by half the present length if the latter
*   is inadequate to store another pair
*
            IF (NVINTI .GT. NDVIN) THEN
               NEWSIZ = NDVIN+NDVIN/2
               CALL RALLOC (PVINIL,NDVIN,NEWSIZ,4)
               CALL RALLOC (PVINVL,NDVIN,NEWSIZ,8)
               NDVIN = NEWSIZ
            ENDIF
*
            DO 2 I = NVINTI,LOC+2,-1
               INDTEI(I) = INDTEI(I-1)
               VALTEI(I) = VALTEI(I-1)
    2       CONTINUE
*
*   Put the new index and value into storage
*
            INDTEI(LOC+1) = INDEX
            VALTEI(LOC+1) = TEGRAL
*
         ENDIF
*
      ELSE
*
*   This branch is executed only once
*
         FIRST = .FALSE.
*
*   Allocate the initial storage for arrays INDTEI and VALTEI;
*   NDVIN stores the array dimension
*
         NDVIN = 2
         CALL ALLOC (PVINIL,NDVIN,4)
         CALL ALLOC (PVINVL,NDVIN,8)
*
*   Compute the integral's value
*
         TEGRAL = VINTI (IA,IB)
*
*   Initialise the integral counter
*
         NVINTI = 1
*
*   Store the integral and its value
*
         INDTEI(1) = INDEX
         VALTEI(1) = TEGRAL
*
      ENDIF
*
      RETURN
      END
