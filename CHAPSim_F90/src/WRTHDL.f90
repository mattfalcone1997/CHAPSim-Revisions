
       SUBROUTINE ERRHDL(msg,rank)
       IMPLICIT NONE
       
       CHARACTER*(*) msg
       INTEGER       rank
       
   11  FORMAT(A,I3.1,5X,A)
       
       write(*,11) '# Error Msg in MYID ', rank, msg
       
       STOP
       
       END SUBROUTINE ERRHDL

          
       SUBROUTINE CHKHDL(msg,rank)
       IMPLICIT NONE
       
       CHARACTER*(*) msg
       INTEGER       rank
       
   11  FORMAT(A,I3.1,5X,A)
       
       write(*,11) '# MYID ', rank, msg
       
       END SUBROUTINE CHKHDL   
       
           
       SUBROUTINE CHKINTHDL(msg,rank,n)
       IMPLICIT NONE
       
       CHARACTER*(*) msg
       INTEGER       rank
       INTEGER       n
       
   11  FORMAT(A,I3.1,5X,A,5X,I9.1)
       
       write(*,11) '# MYID ', rank, msg,n
       
       END SUBROUTINE CHKINTHDL  
       
            
       SUBROUTINE CHKRLHDL(msg,rank,a)
       IMPLICIT NONE
       
       CHARACTER*(*)        msg
       INTEGER              rank
       DOUBLE PRECISION    a
       
   11  FORMAT(A,I3.1,5X,A,5X,ES13.5)
       
       write(*,11) '# MYID ', rank, msg,a
       
       END SUBROUTINE CHKRLHDL  
