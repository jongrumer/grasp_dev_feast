************************************************************************
*
      subroutine items (ncd, ncmin, ncf, record, ierr)
*
* Purpose:
*   Parse a list of levels. 
*
* Input:
*   ncf, record
*
* Output:
*   ierr = 0 - normal 
*          1 - no state
*         -1 - cannot decode
*         -2 - out of range (1,ncf)
*         -3 - cannot decode format like 8-3
*
* Input/Output:
*   ncd - on input 0 then pccmin not allocated;
*             non 0 then allocated
*         on output ==ncmin
*   ncmin - on input, end point of iccmin();
*           on outpu, end poinf of iccmin().
*
* Note:
*   This complicated ncd/ncmin is due to consideration for possible
*   multiple calls of this routine (block version). This structure,
*   however, can be greatly simplified if sufficient memories are
*   been allocated before going into this routine.
*
*   Extracted and modified from rscf92/getold.f  
*   by  Xinghong He                                        Jul 17 1997
************************************************************************
      character record*(*), reci
      CHARACTER*7 FORM, rec*256
      CHARACTER*3 CNUM

      pointer (pccmin,iccmin(*))
      common/def7/pccmin,notused1,notused2
      COMMON/iounit/istdi,istdo,istde

      ncmin_in = ncmin

      if (ncd .eq. 0) then
         ncd = 1
         call alloc (pccmin, ncd, 4)
      endif
*
*   parse record from left to right
*   
      ifirst = 0
      istart = 1
      i = 1
*
*   Remove leading and tailing blanks and obtain the string length
*
      rec = trim ( adjustl (record) )
      length = len_trim (rec)

      i = 1
      do
         reci = rec(i:i)
         if (reci .ge. '0' .and. reci .le. '9') then
            istart = i
            iend = i
         else
            ierr = -1
            return
         endif

         do
            i = i + 1
            if (i .gt. 256) return
            reci = rec(i:i)
            if (reci .ge. '0' .and. reci .le. '9') then
               iend = i
               cycle
            else
               
         

    2 reci = record(i:i)
      if ((reci .ne. ' ') .and. (reci .ne. ',')) then
         istart = i
      else
         i = i+1
         if (i .le. 256) then
            goto 2
         else
            goto 4
         endif
      endif
*
*   .. search for end of string (blank, comma, or dash)
*
    3 reci = record(i:i)
      if ((reci .ne. ' ') .and. (reci .ne. ',') .and.
     :    (reci .ne. '-'))  then
         i = i+1
         if (i .le. 256) goto 3
      endif
*
*     ... read integer
*
      iend = i-1
      isize = iend-istart+1
      call convrt (isize,cnum,lenth)
      form = '(1i'//cnum(1:lenth)//')'
      read (record(istart:iend),form,iostat = ios) level
      if (ios .ne. 0) then
          write(istde,*) 'getold: unable to decode '//
     &                    record(istart:iend)//';'
         ierr = -1
         return
      endif

      if (ifirst .eq. 0) then 
*       .. this is the either the first or an isolated level
         ncmin = ncmin+1
         if (ncmin .gt. ncd) then
             call ralloc (pccmin,ncd,ncmin,4)
             ncd = ncmin
         endif

         if ((level .lt. 1) .or. (level .gt. ncf)) then
            write(istde,*) 'getold: serial numbers must be'
     &,              ' in the range [1,',ncf,'];'
            ierr = -2
            return
         endif

         iccmin(ncmin) = level
         i = i + 1
         if (reci .eq. '-') ifirst = ncmin
         goto 2
      else
*        .. the previous level was the beginning of a range
         level1 = iccmin( ncmin)
         number = level - level1 

         if (number .lt. 0) then
            write(istde,*) level1,'-',level,' not allowed'
            ierr = -3
            return
         endif

         ncmin = ncmin + number
         if (ncmin .gt. ncd) then
             call ralloc (pccmin,ncd,ncmin,4)
             ncd = ncmin
         endif
         do j = 1, number
            iccmin(ifirst+j) =  level1 + j
         enddo
         i = i+1
         ifirst = 0
         goto 2
      endif     

*   at least one level must be requested

    4 if (ncmin .eq. ncmin_in) then
         ierr = 1
         return
      endif
*
*   trim array to exactly the correct size
*
      if (ncmin .ne. ncd) then
         write(istde,*) 'items: ncmin .ne. ncd'
         call ralloc (pccmin,ncd,ncmin,4)
      endif

      ierr = 0
      return
      end
