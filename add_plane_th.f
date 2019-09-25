        Program add_data
        implicit none
        integer        i,j,m,NY,NY_min,NY_max,COUNTER,
     &      TIME_STEP,T_S,T_S_bin,NX,n_pc,NY_TOT,k,NZ
        
        REAL*8        TIME
        LOGICAL       LES 
        parameter     (n_pc = 32,  NY =20, NX=128, NZ=32, 
     &                NY_TOT = NY*n_pc-(n_pc-1),
     &                TIME_STEP=47600,
     &                T_S_bin=1000)

        REAL*8        U1(1:NZ,1:NY_TOT),
     &                U2(1:NZ,1:NY_TOT),
     &                R1(1:NX,1:NY_TOT),
     &                R2(1:NX,1:NY_TOT),
     &                R3(1:NX,1:NY_TOT)
 
        REAL*8        TH1(1:NX,1:NZ)

        CHARACTER*42 filn_th_xy, filn_th_yz, fname_th_xy
        CHARACTER*42 filn_th_xz_1,filn_th_xz_2,filn_th_xz_3


        COUNTER = 0
        open(401,file='add_plane/trans/data_th_xy',status='unknown',
     &           form='unformatted', position='append')
   
        open(402,file='add_plane/trans/data_th_yz',status='unknown',
     &           form='unformatted', position='append')


        open(501,file='add_plane/trans/data_th_xz_1',status='unknown',
     &           form='unformatted', position='append')

        open(502,file='add_plane/trans/data_th_xz_2',status='unknown',
     &           form='unformatted', position='append')

        open(503,file='add_plane/trans/data_th_xz_3',status='unknown',
     &           form='unformatted', position='append')

        open(445,file='add_plane/trans/count_th.time',status='unknown',
     +  form='formatted', position='append')

        open(199,file='add_proc/count_th.dat',status='unknown',
     +  form='formatted', position='append')        
   

        DO T_S=T_S_bin,TIME_STEP,20
         
         COUNTER = COUNTER + 1
 
         fname_th_xy = 'add_proc/data_th_xy_'
     &        //CHAR(MOD(COUNTER,10000000)/1000000+48)
     &        //CHAR(MOD(COUNTER,1000000)/100000+48)
     &        //CHAR(MOD(COUNTER,100000)/10000+48)
     &        //CHAR(MOD(COUNTER,10000)/1000+48)
     &        //CHAR(MOD(COUNTER,1000)/100+48)
     &        //CHAR(MOD(COUNTER,100)/10+48)
     &        //CHAR(MOD(COUNTER,10)+48) //
     &        '.pln' 

         open(603,file=fname_th_xy,status='unknown',
     +  form='unformatted', position='append')

         DO m = 1,n_pc
         
          filn_th_xy = 'plane_data/data_th_xy'
     &         //CHAR(MOD(m,100)/10+48)
     &        //CHAR(MOD(m,10)+48) // '_'
     &        //CHAR(MOD(T_S,10000000)/1000000+48)
     &        //CHAR(MOD(T_S,1000000)/100000+48)
     &        //CHAR(MOD(T_S,100000)/10000+48)
     &        //CHAR(MOD(T_S,10000)/1000+48)
     &        //CHAR(MOD(T_S,1000)/100+48)
     &        //CHAR(MOD(T_S,100)/10+48)
     &        //CHAR(MOD(T_S,10)+48) //
     &        '.pln'

          filn_th_yz = 'plane_data/data_th_yz'
     &         //CHAR(MOD(m,100)/10+48)
     &        //CHAR(MOD(m,10)+48) // '_'
     &        //CHAR(MOD(T_S,10000000)/1000000+48)
     &        //CHAR(MOD(T_S,1000000)/100000+48)
     &        //CHAR(MOD(T_S,100000)/10000+48)
     &        //CHAR(MOD(T_S,10000)/1000+48)
     &        //CHAR(MOD(T_S,1000)/100+48)
     &        //CHAR(MOD(T_S,100)/10+48)
     &        //CHAR(MOD(T_S,10)+48) //
     &        '.pln'  
        
           open(443,file=filn_th_xy,status='old',form='unformatted')

           open(444,file=filn_th_yz,status='old',form='unformatted') 

          NY_min = (m-1)*NY  + 1 - (m-1)
          NY_max = m*NY  - (m-1)
     
          read(443)TIME
          read(443) ((R1(i,j),i=1,NX),j=NY_min,NY_max),
     &              ((R2(i,j),i=1,NX),j=NY_min,NY_max)      
  
         close(443)

          read(444)TIME
          read(444) ((U1(k,j),k=1,NZ),j=NY_min,NY_max),
     &              ((U2(k,j),k=1,NZ),j=NY_min,NY_max)

         close(444)
       


         If (m.eq.1)then
       
          filn_th_xz_1 = 'plane_data/data_th_xz_1'
     &         //CHAR(MOD(m,100)/10+48)
     &        //CHAR(MOD(m,10)+48) // '_'
     &        //CHAR(MOD(T_S,10000000)/1000000+48)
     &        //CHAR(MOD(T_S,1000000)/100000+48)
     &        //CHAR(MOD(T_S,100000)/10000+48)
     &        //CHAR(MOD(T_S,10000)/1000+48)
     &        //CHAR(MOD(T_S,1000)/100+48)
     &        //CHAR(MOD(T_S,100)/10+48)
     &        //CHAR(MOD(T_S,10)+48) //
     &        '.pln'


         open(545,file=filn_th_xz_1,status='old',form='unformatted')

         read(545)TIME
         read(545)((TH1(i,k),i=1,NX),k=1,NZ)

         close(545)

         write(501)TIME
         write(501)((TH1(i,k),i=1,NX),k=1,NZ)

         elseif (m .eq. 3 ) then

           filn_th_xz_2 = 'plane_data/data_th_xz_2'
     &         //CHAR(MOD(m,100)/10+48)
     &        //CHAR(MOD(m,10)+48) // '_'
     &        //CHAR(MOD(T_S,10000000)/1000000+48)
     &        //CHAR(MOD(T_S,1000000)/100000+48)
     &        //CHAR(MOD(T_S,100000)/10000+48)
     &        //CHAR(MOD(T_S,10000)/1000+48)
     &        //CHAR(MOD(T_S,1000)/100+48)
     &        //CHAR(MOD(T_S,100)/10+48)
     &        //CHAR(MOD(T_S,10)+48) //
     &        '.pln'


         open(546,file=filn_th_xz_2,status='old',form='unformatted')

         read(546)TIME
         read(546)((TH1(i,k),i=1,NX),k=1,NZ)

         close(546)

         write(502)TIME
         write(502)((TH1(i,k),i=1,NX),k=1,NZ)  
       

         elseif (m .eq. 4 ) then
         
          filn_th_xz_3 = 'plane_data/data_th_xz_3'
     &         //CHAR(MOD(m,100)/10+48)
     &        //CHAR(MOD(m,10)+48) // '_'
     &        //CHAR(MOD(T_S,10000000)/1000000+48)
     &        //CHAR(MOD(T_S,1000000)/100000+48)
     &        //CHAR(MOD(T_S,100000)/10000+48)
     &        //CHAR(MOD(T_S,10000)/1000+48)
     &        //CHAR(MOD(T_S,1000)/100+48)
     &        //CHAR(MOD(T_S,100)/10+48)
     &        //CHAR(MOD(T_S,10)+48) //
     &        '.pln'


         open(547,file=filn_th_xz_3,status='old',form='unformatted')

         read(547)TIME
         read(547)((TH1(i,k),i=1,NX),k=1,NZ)

         close(547)

         write(503)TIME
         write(503)((TH1(i,k),i=1,NX),k=1,NZ)

         endif





! closing of the processor loop
         ENDDO

          write(199,*)COUNTER, T_S, TIME

          write(401)TIME
          write(401) ((R1(i,j),i=1,NX),j=1,NY_TOT),
     &               ((R2(i,j),i=1,NX),j=1,NY_TOT)

          write(402)TIME
          write(402) ((U1(k,j),k=1,NZ),j=1,NY_TOT),
     &               ((U2(k,j),k=1,NZ),j=1,NY_TOT)

        write(603)TIME
        write(603) ((R1(i,j),i=1,NX),j=1,NY_TOT),
     &             ((R2(i,j),i=1,NX),j=1,NY_TOT)
 
        close(603)  
           

        IF (MOD(COUNTER,10).EQ.0) THEN
           WRITE(6,*) COUNTER, T_s, TIME   
        END IF         
        
! closing of the time step  loop
        ENDDO

        write(445,*) counter
        close(445)
        stop
        end



