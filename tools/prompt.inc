c
c Prompt user for input. Associated with USE_LIST_OPT. All RD_LIST_OPT options
c use thus same section of code.
c
c Subroutine may be better, but requires some manipulation.
c
c Created 01-May-2024
c
      if(len_trim(description) .gt. 48)write(*,"(1x,a)")trim(description)
 300  if(len_trim(description) .gt. 48)then
        if(len_trim(default).ne.0)then
          write(*,"(t50,2a,' [',a,']=>')",advance="no")
     *       ' :: ',var_name,trim(default)
         else
            write(*,"(t50,2a,' =>')",advance="no")' :: ',var_name
         endif
       else
        if(len_trim(default).ne.0)then
          write(*,"(1x,a,t50,2a,' [',a,']=>')",advance="no")
     *       trim(description),' :: ',var_name,trim(default)
         else
            write(*,"(1x,a,t50,2a,' =>')",advance="no")trim(description),' :: ',var_name
         endif
       end if
      read(*,"(a)")answer

