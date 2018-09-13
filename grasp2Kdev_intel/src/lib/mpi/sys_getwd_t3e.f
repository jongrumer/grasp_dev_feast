      subroutine sys_getwd (dir,len_cwd) 
      CHARACTER (len=*), intent(out):: dir
      INTEGER GETCWD, len_cwd

      len_cwd = GETCWD(dir)

      if (len_cwd.le. 0) then
         print *, 'Can''t get the current directory: exit!'
         call exit(1);
      endif

      end 
        
