!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
program rcsfzerofirst
!
!  Written by Per Jonsson, Malmo University, Sweden
!
!  Updated January 2012
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

implicit none

integer             :: i, j, k, ndup, nzero(20), nfirst(20)
character(len=180) :: csf11(2000000), csf12(2000000), csf13(2000000)
character(len=180) :: csf21, csf22, csf23
character(len=1000) :: string


write(*,*)
write(*,*)  ' RCSFZEROFIRST'
write(*,*)  ' This program takes a list of CSFs and partitions each'
write(*,*)  ' symmetry block into a zero- and a first-order CSF space'
write(*,*)  ' The CSFs that define the zero-order spaces are read from'
write(*,*)  ' from a zero-order list'
write(*,*)  ' Input files: list with CSFs to be partitioned'
write(*,*)  '              list with CSFs defining the zero-order space'
write(*,*)  ' Output file: rcsf.out'
write(*,*)
write(*,*)  ' Current version allows only blocks in the zero-order space'
write(*,*)  ' less than 2 000 000 CSFs'
write(*,*)

! Get the names of the lists and open the files

write(*,*) ' Give the full name of the list that contains the zero-order spaces'
read(*,'(a)') string
open(19, file = trim(string), status = 'old')

write(*,*) ' Give the full name of the list that should be partitioned'
read(*,'(a)') string
open(20, file = trim(string), status = 'old')

! Open the output file

open(21, file = 'rcsf.out', status = 'unknown')

! Read header of list1

do i = 1,5
  read(19,'(a)') string
end do

! Read header of list2 and write to output file

do i = 1,5
  read(20,'(a)') string
  write(21,'(a)') trim(string)
end do

! Loop over blocks (both list should have the same number blocks)

i = 1

do

! Loop over CSFs in a current block

  j = 1

  do

! Read CSFs in list1 and write to rcsl.out

    read(19,'(a)',end = 99) string

    if (string(1:2) .eq. ' *') then
      exit
    else
      csf11(j) = string
      read(19,'(a)') csf12(j)
      read(19,'(a)') csf13(j)
      write(21,'(a)') trim(csf11(j))
      write(21,'(a)') trim(csf12(j))
      write(21,'(a)') trim(csf13(j))
    end if

    j = j + 1
    if (j.gt.2000000) then
       write(*,*) 'In rscfzerofirst'
       write(*,*) 'Increase size of csf11, csf12, csf13'
       write(*,*) 'and recompile'
       STOP
    end if

  end do

99 continue

  nzero(i) = j - 1

  nfirst(i) = 0

  do

! Read CSFs in list2 and see if it is in list1
! If not write to rcsl.out

    read(20,'(a)',end = 999) string

    if (string(1:2) .eq. ' *') then
      exit
    else
      csf21 = string
      read(20,'(a)') csf22
      read(20,'(a)') csf23
    end if

    ndup = 0
    do k = 1,nzero(i)
      if ((csf11(k) == csf21).and.(csf12(k) == csf22).and.(csf13(k) == csf23)) then
        ndup = 1
        exit
      end if
    end do    
    if (ndup == 0) then
      write(21,'(a)') trim(csf21)
      write(21,'(a)') trim(csf22)
      write(21,'(a)') trim(csf23)
      nfirst(i) = nfirst(i) + 1
    end if
  end do

! write * to mark end of block

  write(21,'(a)') ' *'

  i = i + 1

end do

999 continue

! Write out information about the blocks and the zero-order space

write(*,*)
write(*,*) 'Block,   CSFs in zero-order space,     Total number of CSFs'

do j = 1,i
  write(*,200) j,nzero(j),nzero(j) + nfirst(j)
end do
write(*,*) 
write(*,*) 'Partitioned list written to rcsf.out'

200 format(i4,i20,i25)

end program rcsfzerofirst
