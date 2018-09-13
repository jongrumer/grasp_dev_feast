************************************************************************
*                                                                      *
      SUBROUTINE GENMCP(NAME,IC,NTESTG,INPCI)
*                                                                      *
*   Entry routine for GENMCP. Controls the computation of the          *
*   one-particle coupling coefficients.                                *
*                                                                      *
*   Written by Per Jonsson                                 June 1996   *
*                                                                      *
************************************************************************
*
      CHARACTER*24 NAME
*
      PRINT *
      PRINT *, 'GENMCP: Execution begins for ',NAME
*
*   Open, check, load data from, and close the  csl  file
*
      CALL SETCSLA(NAME, ncore_not_used)
*
*   Set up the table of logarithms of factorials for use by
*   angular modules
*
      CALL FACTT
*
*   Proceed with the generation of MCP coefficients
*
      CALL MCPOUT(NAME,IC,NTESTG,INPCI)
      CALL MCPIN(NAME,IC,NTESTG,INPCI)
*
*   Print completion message
*
      PRINT *
*
      RETURN
      END
