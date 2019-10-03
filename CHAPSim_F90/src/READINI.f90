       SUBROUTINE READINI
       use init_info
       use mpi_info
       use mesh_info
       IMPLICIT NONE
       
       CHARACTER(128) :: sect         
       CHARACTER(128) :: skey          
       INTEGER(4) :: IOS=0
       INTEGER(4) :: INI=13           
       INTEGER(4) :: lens    
       
       OPEN(INI,FILE='readdata.ini', status='old', iostat=ios)
       if(ios .NE. 0)  &
       call ERRHDL(' File readdata.ini cannot be found.',MYID)
      
       read(ini,*) sect
       do while(sect(1:1).EQ.';' .or. sect(1:1).EQ.'#'.or. &
          sect(1:1).EQ.' ')
          read(ini,*) sect   
       end do
       lens = len_trim(sect)
       call CHKHDL(' Read in '//sect(1:lens),MYID)
       
  
       if(sect(1:lens) /= '[geometry]')  &
       call ERRHDL(' Reading fails: '//sect(1:lens),MYID)

       read(ini,*) skey, HX
       call CHKRLHDL  ('   HX=          ',MYID,HX)
       read(ini,*) skey, HZ
       call CHKRLHDL  ('   HZ=          ',MYID,HZ)
       read(ini,*) skey, iswitch

       IF(iswitch==1) call CHKHDL ('   iswitch=    Channel Flow',MYID) 
       IF(iswitch==2) call CHKHDL ('   iswitch=    Pipe Flow',MYID) 
       read(ini,*) skey, STR2
       call CHKRLHDL  ('   STR2=        ',MYID,STR2)
       read(ini,*) skey, ISTR2
       call CHKINTHDL ('   ISTR2=       ',MYID,ISTR2)
         
       read(ini,*) sect
       lens = len_trim(sect)
       call CHKHDL(' Read in '//sect(1:lens),MYID)  
       if(sect(1:lens) /= '[mesh]')       &
       call ERRHDL(' Reading fails: '//sect(1:lens),MYID)     
       
       read(ini,*) skey, NCL1,NCL1_io
       call CHKINTHDL ('   NCL1= (MAIN) ',MYID,NCL1_io )
       call CHKINTHDL ('   NCL1= (TurGe)',MYID,NCL1 )
       read(ini,*) skey, NCL2
       call CHKINTHDL ('   NCL2=        ',MYID,NCL2)
       read(ini,*) skey, NCL3
       call CHKINTHDL ('   NCL3=        ',MYID,NCL3)
       
       IF(NCL1_io .GT. 2) THEN
          IOFLOWflg =  .true. 
       ELSE
          IOFLOWflg =  .false. 
       END IF
       
       read(ini,*) sect
       lens = len_trim(sect)
       call CHKHDL(' Read in '//sect(1:lens),MYID)  
       if(sect(1:lens) /= '[Pboundary]')       &
       call ERRHDL(' Reading fails: '//sect(1:lens),MYID)     
       
       read(ini,*) skey, BCX(1)
       call CHKINTHDL  ('   BCX(1)=      ',MYID,BCX(1))
       read(ini,*) skey, BCX(2)
       call CHKINTHDL  ('   BCX(2)=      ',MYID,BCX(2))
       if( BCX(1) == 3 ) BCX(2) = 3
       if( BCX(2) == 3 ) BCX(1) = 3
       
       read(ini,*) skey, BCX_io(1)
       call CHKINTHDL  ('   BCX_io(1)=   ',MYID,BCX_io(1))
       read(ini,*) skey, BCX_io(2)
       call CHKINTHDL  ('   BCX_io(2)=   ',MYID,BCX_io(2))
       if( BCX_io(1) == 3 ) BCX_io(2) = 3
       if( BCX_io(2) == 3 ) BCX_io(1) = 3
       
       
       read(ini,*) skey, BCZ(1)
       call CHKINTHDL  ('   BCZ(1)=      ',MYID,BCZ(1))
       read(ini,*) skey, BCZ(2)
       call CHKINTHDL  ('   BCZ(2)=      ',MYID,BCZ(2))
       if( BCZ(1) == 3 ) BCZ(2) = 3
       if( BCZ(2) == 3 ) BCZ(1) = 3
       
       IF(BCZ(1)/=3) &
       call ERRHDL('Z MUST BE PERIODIC',MYID)  

         
       read(ini,*) sect
       lens = len_trim(sect)
       call CHKHDL(' Read in '//sect(1:lens),MYID)  
       if(sect(1:lens) /= '[numerical]')       &
       call ERRHDL(' Reading fails: '//sect(1:lens),MYID)     
       
       read(ini,*) skey, DT
       call CHKRLHDL  ('   DT=       ',MYID,DT)
       read(ini,*) skey, MULTIM
       call CHKINTHDL ('   MULTIM=   ',MYID,MULTIM)
       read(ini,*) skey, CFLGV
       call CHKRLHDL  ('   CFLGV=    ',MYID,CFLGV)
       read(ini,*) skey, NREAD, NREAD_io
       IF(NREAD==0) call CHKHDL ('   NREAD=      from random flow field',MYID)
       IF(NREAD==1) call CHKHDL ('   NREAD=      interpolation from a coarse mesh',MYID)
       IF(NREAD==2) call CHKHDL ('   NREAD=      restart from last step',MYID)
       IF(NREAD_io==0) call CHKHDL ('   NREAD_io=   from random flow field',MYID)
       IF(NREAD_io==1) call CHKHDL ('   NREAD_io=   interpolation from a coarse mesh',MYID)
       IF(NREAD_io==2) call CHKHDL ('   NREAD_io=   restart from last step',MYID)
       
       read(ini,*) skey, TRST, TRST_io
       call CHKRLHDL  ('   TRST=     ',MYID,TRST)
       call CHKRLHDL  ('   TRST_io=  ',MYID,TRST)
       read(ini,*) skey, TSTOP 
       call CHKRLHDL  ('   TSTOP=    ',MYID,TSTOP)
       read(ini,*) skey, tbody
       call CHKRLHDL  ('   tbody=    ',MYID,tbody)
       read(ini,*) skey, TLGRE
       call CHKRLHDL  ('   TLGRE=    ',MYID,TLGRE)
       read(ini,*) skey, REINI
       call CHKRLHDL  ('   REINI=    ',MYID,REINI)
       
          
       read(ini,*) sect
       lens = len_trim(sect)
       call CHKHDL(' Read in '//sect(1:lens),MYID)  
       if(sect(1:lens) /= '[fluid]') &
       call ERRHDL(' Reading fails: '//sect(1:lens),MYID) 
       
       read(ini,*) skey, NFLOW 
       if(NFLOW==1)call CHKHDL ('   NFLOW=      X Streamwise Flow',MYID)
       if(NFLOW==2)call CHKHDL ('   NFLOW=      Y Streamwise Flow',MYID)
       if(NFLOW==3)call CHKHDL ('   NFLOW=      Z Streamwise Flow',MYID)
       read(ini,*) skey, REN
       call CHKRLHDL  ('   REN=      ',MYID,REN)
       read(ini,*) skey, VPER
       call CHKRLHDL  ('   VPER=     ',MYID,VPER)
       read(ini,*) skey, SVPER
       call CHKRLHDL  ('   SVPER=    ',MYID,VPER)
       read(ini,*) skey, FLOWTP
       if(FLOWTP==1) call CHKHDL ('   FLOWTP=     Constant mass flux driven',MYID)
       if(FLOWTP==2) call CHKHDL ('   FLOWTP=     Constant pressure gradient driven',MYID)
       read(ini,*) skey, CFGV
       call CHKRLHDL  ('   CFGV=    ',MYID,CFGV)
       read(ini,*) skey, RATEM1
       call CHKRLHDL  ('   RATEM1=   ',MYID,RATEM1)

       read(ini,*) sect
       lens = len_trim(sect)
       call CHKHDL(' Read in '//sect(1:lens),MYID)  
       if(sect(1:lens) /= '[statistics]')  &
       call ERRHDL(' Reading fails: '//sect(1:lens),MYID) 
       
       read(ini,*) skey, TSTAV1
       call CHKRLHDL  ('   TSTAV1=   ',MYID,TSTAV1)
       read(ini,*) skey, TTECCK
       call CHKRLHDL  ('   TTECCK=   ',MYID,TTECCK)
       read(ini,*) skey, TSAVE1
       call CHKRLHDL  ('   TSAVE1=   ',MYID,TSAVE1)
       read(ini,*) skey, TSCN
       call CHKRLHDL  ('   TSCN=     ',MYID,TSCN)
       read(ini,*) skey, N_START
       call CHKINTHDL ('   N_START=  ',MYID,N_START)
       
       close(ini)

       NND1=NCL1+1
       NND2=NCL2+1
       NND3=NCL3+1
       
       NND1_io=NCL1_io+1                   
       
       RETURN
             
       END SUBROUTINE READINI
       
       
       
       
       
