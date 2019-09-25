


 
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE SAVE_STATS_CHAN(FINAL)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi' 

      CHARACTER*35 FNAME
      CHARACTER*35 FNAME_TH
      LOGICAL FINAL
      integer i,j,k,n
      real*8 uc, ubulk
      real*16 vblf
      real*8 PV(0:NY+1)
!      real*8 GYF_TOT(1:NY*NPROCS),UME_TOT(1:NY*NPROCS),
!     +       WME_TOT(1:NY*NPROCS),VME_TOT(1:NY*NPROCS),TH_TOT(1:NY*NPRO
!     + CS), SAL_TOT(1:NY*NPROCS),urms_tot(1:NY*NPROCS),
!     +       vrms_tot(1:NY*NPROCS),wrms_tot(1:NY*NPROCS),
!     +       vbl(1:NY*NPROCS),shear_t(1:NY*NPROCS),
!     +       dthdy_t(1:NY*NPROCS),dsdy_t(1:NY*NPROCS),
!     +       thrms_t(1:NY*NPROCS),srms_t(1:NY*NPROCS),
!     +       tke_t(1:NY*NPROCS),prod_t(1:NY*NPROCS),
!     +       buo(1:NY*NPROCS),bth(1:NY*NPROCS),bs(1:NY*NPROCS),
!     +       diss_t(1:NY*NPROCS),uv_t(1:NY*NPROCS),
!     +       wv_t(1:NY*NPROCS),uw_t(1:NY*NPROCS)

      if (rank .eq.0)then
      WRITE(6,*) 'Saving flow statistics.'
      endif

      IF (USE_MPI) THEN
          CALL GHOST_CHAN_MPI
          CALL  GHOST_ZERO_CHAN_MPI
      ENDIF

        I = RANK+1

           FNAME ='results/mean_'
     &        //CHAR(MOD(I,100)/10+48)
     &        //CHAR(MOD(I,10)+48) // '.txt'

           FNAME_TH ='results/mean_th_'
     &        //CHAR(MOD(I,100)/10+48)
     &        //CHAR(MOD(I,10)+48) // '.txt'
        

c       write(400,*) time_step, CTH(0,0,1,1),  10000
c       write(410,*) time_step, CTH(0,0,NY,1), 20000
c       write(500,*) time_step, CU1(0,0,0),CU1(0,0,1)

C Apply Boundary conditions to velocity field
      IF (USE_MPI) THEN
        CALL APPLY_BC_VEL_MPI
      ELSE
        IF ((U_BC_YMIN.EQ.3).or.(W_BC_YMIN.EQ.3)) THEN
          CALL APPLY_BC_NWM_LOWER
        ELSE
          CALL APPLY_BC_VEL_LOWER
        ENDIF

        IF ((U_BC_YMAX.EQ.3).or.(W_BC_YMAX.EQ.3)) then
          CALL APPLY_BC_NWM_UPPER
        ELSE
          CALL APPLY_BC_VEL_UPPER
        ENDIF
      END IF


! We are in the middle of a run, compile statistics
! First get the number of samples taken so far
      IF (RANK .eq. 0)THEN
      write(*,*) 'TIME, DELTA_T: ',TIME, DELTA_T
      ENDIF
! Compute and write out bulk velocity
! Integrat the instantaneous mean profile numerically at GY points
      UBULK=0.
      do j=2,NY
        UBULK=UBULK+0.5*(dble(CU1(0,0,j)+CU1(0,0,j-1)))*DY(j)
      end do
      UBULK=UBULK/LY      
! Write out UBULK
      if (rank .eq.0)then
      write(*,*) 'UBULK: ',UBULK
      endif

! Save CUi
      do k=0,TNKZ
        do i=0,NKX
          do j=0,NY+1
            CR1(i,k,j)=CU1(i,k,j)
            CR2(i,k,j)=CU2(i,k,j)
            CR3(i,k,j)=CU3(i,k,j)
            CFTH(i,k,j,1)=CP(i,k,j)
          end do
        end do
      end do
    
c      write(400,*) time_step, CU1(0,0,0),  10000
c      write(410,*) time_step, CU1(0,0,NY), 20000
 
! Convert to physical space
c      call fft_xz_to_physical(CU1(0,0,0),U1(0,0,0),0,NY+1)
c      call fft_xz_to_physical(CU2(0,0,0),U2(0,0,0),0,NY+1)
c      call fft_xz_to_physical(CU3(0,0,0),U3(0,0,0),0,NY+1)

      CALL FFT_XZ_TO_PHYSICAL(CU1(0,0,-1),U1(0,0,-1),0,NY+1)
      CALL FFT_XZ_TO_PHYSICAL(CU2(0,0,-1),U2(0,0,-1),0,NY+1)
      CALL FFT_XZ_TO_PHYSICAL(CU3(0,0,-1),U3(0,0,-1),0,NY+1)

      CALL FFT_XZ_TO_PHYSICAL(CU1(0,0,NY+1),U1(0,0,NY+1),0,1)
      CALL FFT_XZ_TO_PHYSICAL(CU2(0,0,NY+1),U2(0,0,NY+1),0,1)
      CALL FFT_XZ_TO_PHYSICAL(CU3(0,0,NY+1),U3(0,0,NY+1),0,1)

      CALL fft_xz_to_physical(CP(0,0,0),P(0,0,0),0,NY+1)

c       call tkebudget_chan 

c       open(54,file='field_data_xz_49.txt')
c       open(55,file='field_data_xz_98.txt')

c      if (MOD(TIME_STEP,50).EQ.0) then
c          write(54,*) TIME_STEP, TIME, dble(CR1(0,0,NY))
c          write(55,*) TIME_STEP, TIME, dble(CR1(0,0,NY))
c      endif


! Get the turbulent kinetic energy at each level 
      do j=0,NY+1
        urms(j)=0.
        vrms(j)=0.
        wrms(j)=0.
      do k=0,NZM
      do i=0,NXM 
        urms(j)=urms(j)+(abs(U1(i,k,j)-dble(CR1(0,0,j))))**2.
!        urms(j)=urms(j)+(abs(U1(i,k,j)))**2.
        vrms(j)=vrms(j)+(abs(U2(i,k,j)-dble(CR2(0,0,j))))**2.
!        vrms(j)=vrms(j)+(abs(U2(i,k,j)))**2.
        wrms(j)=wrms(j)+(abs(U3(i,k,j)-dble(CR3(0,0,j))))**2.
!         wrms(j)=wrms(j)+(abs(U3(i,k,j)))**2.

c        if ((MOD(TIME_STEP,50).EQ.0).and. (j .eq. NY/3)) then
c          write(54,*) U1(i,k,j)-dble(CR1(0,0,j)),
c     :        U2(i,k,j)-dble(CR2(0,0,j)),
c     :        U3(i,k,j)-dble(CR3(0,0,j)),U1(i,k,j),U2(i,k,j),U3(i,k,j)
c        elseif ((MOD(TIME_STEP,50).EQ.0).and.(j .eq. NY/2)) then
c          write(55,*) U1(i,k,j)-dble(CR1(0,0,j)),
c     :        U2(i,k,j)-dble(CR2(0,0,j)),
c     :        U3(i,k,j)-dble(CR3(0,0,j)),U1(i,j,k),U2(i,j,k),U3(i,j,k)
c        endif
C         if(k==10) then
         urms2d(i,j)=urms2d(i,j)+(abs(U1(i,k,j)-dble(CR1(0,0,j))))**2.
         vrms2d(i,j)=vrms2d(i,j)+(abs(U2(i,k,j)-dble(CR2(0,0,j))))**2.
         wrms2d(i,j)=wrms2d(i,j)+(abs(U3(i,k,j)-dble(CR3(0,0,j))))**2.
C         endif
      end do
      end do
        urms(j)=sqrt(urms(j)/(float(NZ)*float(NX)))
!        urms1(j)=sqrt(urms(j)/(float(NZ)*float(NX)))
        vrms(j)=sqrt(vrms(j)/(float(NZ)*float(NX)))
!        vrms1(j)=sqrt(vrms(j)/(float(NZ)*float(NX)))
        wrms(j)=sqrt(wrms(j)/(float(NZ)*float(NX)))
!        wrms1(j)=sqrt(wrms(j)/(float(NZ)*float(NX)))
      end do 
! Get the bulk rms value
      do i=1,NX
      do j=0,NY+1
      write(668,*) wrms2d(i,j)
      urms2d(i,j)=sqrt(urms2d(i,j))!/float(NZ))
      vrms2d(i,j)=sqrt(vrms2d(i,j))!/float(NZ))
      wrms2d(i,j)=sqrt(wrms2d(i,j))!/float(NZ))
      enddo
      enddo
      urms_b=0.
      do j=2,NY
        urms_b=urms_b+0.5*(urms(j)+urms(j-1))*DY(j)
      end do
      urms_b=urms_b/LY      
      vrms_b=0.
      do j=1,NY
        vrms_b=vrms_b+vrms(j)*DY(j)
      end do
      vrms_b=vrms_b/LY
      wrms_b=0.
      do j=2,NY
        wrms_b=wrms_b+0.5*(wrms(j)+wrms(j-1))*DY(j)
      end do
      wrms_b=wrms_b/LY
! Compute the Reynolds stress and mean velocity gradient
      do j=1,NY
        uv(j)=0. 
        uw(j)=0.
        wv(j)=0.
        pv(j)=0.
      do k=0,NZM
      do i=0,NXM
        uv(j)=uv(j)+(U1(i,k,j)-dble(CR1(0,0,j)))
     +    *(0.5*(U2(i,k,j)+U2(i,k,j+1))
     &    -0.5*dble(CR2(0,0,j)+CR2(0,0,j+1)))
        wv(j)=wv(j)+(U3(i,k,j)-dble(CR3(0,0,j)))
     +    *(0.5*(U2(i,k,j)+U2(i,k,j+1))
     &      -0.5*dble(CR2(0,0,j)+CR2(0,0,j+1)))
        uw(j)=uw(j)+(U1(i,k,j)-dble(CR1(0,0,j)))
     +    *(U3(i,k,j)-dble(CR3(0,0,j)))
        pv(j)=pv(j)+(P(i,k,j)-dble(CFTH(0,0,j,1)))
     +    *(0.5*(U2(i,k,j)+U2(i,k,j+1))
     &    -0.5*dble(CR2(0,0,j)+CR2(0,0,j+1)))
      end do
      end do
        uv(j)=uv(j)/(float(NZ)*float(NX))
        uw(j)=uw(j)/(float(NZ)*float(NX))
        wv(j)=wv(j)/(float(NZ)*float(NX))
        pv(j)=pv(j)/(float(NZ)*float(NX))
      end do
              
! Get the y-derivative of the mean velocity at GYF points
!      do j=1,NY
!        dudy(j)=dble(CR1(0,0,j+1)-CR1(0,0,j-1))/(2.*DYF(j))
!        dwdy(j)=dble(CR3(0,0,j+1)-CR3(0,0,j-1))/(2.*DYF(j))
!      end do
! Get the y-derivative of the mean velocity at GY points
      do j=1,NY
        dudy(j)=dble(CR1(0,0,j)-CR1(0,0,j-1))/(GYF(j)-GYF(j-1))
        dwdy(j)=dble(CR3(0,0,j)-CR3(0,0,j-1))/(GYF(j)-GYF(j-1))
      end do
 
c      write(400,*) time_step0, CR1(0,0,0),   10000
c      write(410,*) time_step, CR1(0,0,NY-1), 20000
c       write(500,*) time_step, CR1(0,0,0),CR1(0,0,1)

! Calculate the mean square shear
      do j=1,NY
        shear(j)=0.d0
        do k=0,NZM
          do i=0,NXM
            shear(j)=shear(j)
     &            +((U1(i,k,j+1)-U1(i,k,j-1))/(2.d0*DYF(j)))**2.d0
     &            +((U3(i,k,j+1)-U3(i,k,j-1))/(2.d0*DYF(j)))**2.d0
          end do
        end do
        shear(j)=shear(j)/dble(NX*NZ)
      end do

! Write out the bulk rms velocity
      IF (RANK .eq. 0) THEN 
       write(*,*) '<U_rms>: ',urms_b
       write(*,*) '<V_rms>: ',vrms_b
       write(*,*) '<W_rms>: ',wrms_b
      ENDIF

! Get the rms vorticity
! First, get the x-component in fourier space
      do j=1,NY
      do k=0,TNKZ
      do i=0,NKX
        CS1(i,k,j)=(CR3(i,k,j+1)-CR3(i,k,j-1))/(2.d0*DYF(j))
     &            -CIKZ(K)*0.5d0*(CR2(i,k,j+1)+CR2(i,k,j))
      end do
      end do
      end do
! Convert to physical space
      call fft_xz_to_physical(CS1,S1,0,NY+1)
! Get the rms value
      do j=1,NY
      omega_x(j)=0.d0
      do k=1,NZM
      do i=1,NXM
        omega_x(j)=omega_x(j)+S1(i,k,j)**2.d0
      end do
      end do
      omega_x(j)=sqrt(omega_x(j)/(dble(NX-1)*dble(NZ-1)))
      end do

! Now, get the y-component in fourier space
      do j=1,NY
      do k=0,TNKZ
      do i=0,NKX
        CS1(i,k,j)=CIKZ(k)*CR1(i,k,j)-CIKX(i)*CR3(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xz_to_physical(CS1,S1,0,NY+1)
! Get the rms value
      do j=1,NY
      omega_y(j)=0.d0
      do k=0,NZM
      do i=0,NXM
        omega_y(j)=omega_y(j)+S1(i,k,j)**2.d0
      end do
      end do
      omega_y(j)=sqrt(omega_y(j)/(dble(NX)*dble(NZ)))
      end do

! Now, get the z-component in fourier space
      do j=1,NY
      do k=0,TNKZ
      do i=0,NKX
        CS1(i,k,j)=CIKX(i)*0.5d0*(CR2(i,k,j+1)+CR2(i,k,j))
     &             -(CR1(i,k,j+1)-CR1(i,k,j-1))/(2.d0*DYF(j))
      end do
      end do
!        CS1(0,0,j)=CS1(0,0,j)-dudy
      end do
! Convert to physical space
      call fft_xz_to_physical(CS1,S1,0,NY+1)
! Get the rms value
!      write(555,*) CS1(0,0,1) 
      do j=1,NY
      omega_z(j)=0.d0
      do k=0,NZM
      do i=0,NXM
        omega_z(j)=omega_z(j)+(S1(i,k,j))**2.d0
      end do
      end do
      omega_z(j)=sqrt(omega_z(j)/(dble(NX)*dble(NZ)))
      end do


      if (MOVIE) then
!        call save_plane_vel
      endif

      

! Write out the mean statistics at each time
      open(40,file=FNAME,form='formatted',status='unknown',
     + position='append')

      write(40,*) TIME_STEP,TIME,DELTA_T,UBULK
      do j=1,NY
        write(40,411) j,GYF(J),dble(CR1(0,0,j))
     +      ,0.5*dble((CR2(0,0,j-1)+CR2(0,0,j)))
     +      ,dble(CR3(0,0,j)),urms(j),vrms(j),wrms(j)
     +      ,uv(j),uw(j),wv(j),dudy(j),dwdy(j),dble(CFTH(0,0,j,1))
     +      ,shear(j),omega_x(j),omega_y(j),omega_z(j),pv(j)
      end do
      close(40)
411   format(I3,' ',19(F20.9,' '))
c      write(400,*) time_step, CU2(0,0,1),  10000
c      write(410,*) time_step, CU2(0,0,NY), 20000

! If we are using a near wall model, write out prescribed wall stress
      if ((U_BC_YMIN.eq.3).or.(W_BC_YMIN.eq.3)
     & .or.(U_BC_YMAX.eq.3).or.(W_BC_YMAX.eq.3)) then
        open(50,file='mean_nwm.txt',form='formatted',status='unknown')
        write(50,*) TIME_STEP,TIME,DELTA_T,UBULK
        write(50,501) SUM(U_BC_LOWER_NWM(0:NXM,0:NZM))/dble(NX*NZ),
     &      SUM(W_BC_LOWER_NWM(0:NXM,0:NZM))/dble(NX*NZ),
     &      SUM(U_BC_UPPER_NWM(0:NXM,0:NZM))/dble(NX*NZ),
     &      SUM(W_BC_UPPER_NWM(0:NXM,0:NZM))/dble(NX*NZ),
     &      tauwall_mean

      end if
501   format(5(F25.9,' '))

 
C      call tkebudget_chan_1   

C      Call netcdf
C       call   NETCDF_WRITE
      thu2d(:,:,:)=0.0d0
      thv2d(:,:,:)=0.0d0

! Do over the number of passive scalars
      do n=1,N_TH

! Save CTH
      do k=0,TNKZ
        do i=0,NKX
          do j=0,NY+1
            CRTH(i,k,j,n)=CTH(i,k,j,n)
          end do
        end do
      end do

! Compute the scalar gradient and store in CRi
      do j=1,NY
        do k=0,TNKZ
          do i=0,NKX
! Store gradients of TH(:,:,:,n) (if it is used) in CRi
          CF1(i,k,j)=CIKX(i)*CTH(i,k,j,n)
          CF2(i,k,j)=(CTH(i,k,j+1,n)-CTH(i,k,j-1,n))/(GYF(j+1)-GYF(j-1))
          CF3(i,k,j)=CIKZ(k)*CTH(i,k,j,n)
          end do
        end do
      end do
! Convert gradients to physical space
      CALL FFT_XZ_TO_PHYSICAL(CF1,F1,0,NY+1)
      CALL FFT_XZ_TO_PHYSICAL(CF2,F2,0,NY+1)
      CALL FFT_XZ_TO_PHYSICAL(CF3,F3,0,NY+1)

! Convert to physical space
      call fft_xz_to_physical(CTH(0,0,0,n),TH(0,0,0,n),0,NY+1)


c      open(56,file='field_data_th_xz_49.txt')
c      open(57,file='field_data_th_xz_98.txt')

c      if (MOD(TIME_STEP,50).EQ.0) then
c        write(56,*)TIME_STEP, TIME
c        write(57,*)TIME_STEP, TIME
c      endif

      do j=1,NY
        thrms(j,n)=0.
      do k=0,NZM
      do i=0,NXM
        thrms(j,n)=thrms(j,n)+(abs(TH(i,k,j,n)-dble(CRTH(0,0,j,n))))**2.

c        if ((MOD(TIME_STEP,50).EQ.0).and.(j .eq. NY/3)) then
c          write(56,*) TH(i,k,j,n)-dble(CRTH(0,0,j,n)), TH(i,k,j,n)
c        elseif ((MOD(TIME_STEP,50).EQ.0).and.(j .eq. NY/2)) then
c          write(57,*) TH(i,k,j,n)-dble(CRTH(0,0,j,n)), TH(i,k,j,n)
c        endif

      end do
      end do
        thrms(j,n)=sqrt(thrms(j,n)/(float(NZ)*float(NX)))
      end do
! Compute the Reynolds stress and mean velocity gradient
      do j=1,NY
        thu(j,n)=0.
        thv(j,n)=0.
      do k=0,NZM
      do i=0,NXM
       thu(j,n)=thu(j,n)+(TH(i,k,j,n)-dble(CRTH(0,0,j,n)))
     +    *(U1(i,k,j) -dble(CR1(0,0,j)))
       thv(j,n)=thv(j,n)+(TH(i,k,j,n)-dble(CRTH(0,0,j,n)))
     +    *(0.5*(U2(i,k,j)+U2(i,k,j+1))
     &      -0.5*dble(CR2(0,0,j)+CR2(0,0,j+1)))
       thu2d(i,j,n)=thu2d(i,j,n)+(TH(i,k,j,n)-dble(CRTH(0,0,j,n)))
     +    *(U1(i,k,j) -dble(CR1(0,0,j)))
       thv2d(i,j,n)=thv2d(i,j,n)+(TH(i,k,j,n)-dble(CRTH(0,0,j,n)))
     +    *(0.5*(U2(i,k,j)+U2(i,k,j+1))
     +      -0.5*dble(CR2(0,0,j)+CR2(0,0,j+1)))
      end do
      end do
      thu(j,n)=thu(j,n)/(float(NZ)*float(NX))
      thv(j,n)=thv(j,n)/(float(NZ)*float(NX))
      if(n>1)then
       thus(j)=thu(j,n)
       thvs(j)=thv(j,n)
      else
       thut(j)=thu(j,n)
       thvt(j)=thv(j,n)
      endif
      end do
      
      do j=1,NY
        do i=0,NXM
         thu2d(i,j,n)=thu2d(i,j,n)/float(NZ)
         thv2d(i,j,n)=thv2d(i,j,n)/float(NZ)
        enddo
      enddo
! Get the y-derivative of the mean scalar at GYF points
      do j=1,NY
        dthdy(j,n)=dble(CRTH(0,0,j+1,n)-CRTH(0,0,j,n))
     &           /(GYF(j+1)-GYF(j))
        if(n>1)then
        dthdys(j)=dthdy(j,n)
        else
        dthdyt(j)=dthdy(j,n)
        endif
      end do

c       write(500,*) time_step, CRTH(0,0,0,1),CRTH(0,0,1,1)

! Compute the potential energy dissipation, grad(TH) \cdot grad(TH)
      do j=1,NY
        pe_diss(j,n)=0.d0
        do k=0,NZM
          do i=0,NXM
            pe_diss(j,n)=pe_diss(j,n)
     &          +F1(i,k,j)**2.d0+F2(i,k,j)**2.d0+F3(i,k,j)**2.d0
          end do
        end do
        pe_diss(j,n)=pe_diss(j,n)/dble(NX*NZ)
      end do

      end do !(loop for N_TH)

       call tkebudget_chan_1
!!!!!!!!!!! Combine all statistics !!!!!!!!!!!!!

       call MPI_PROF_ALL_STAT
     
!      call MPI_INST_PROF(dble(CR1(0,0,1:NY)),UME_TOT)
!      call MPI_INST_PROF(dble(CR2(0,0,1:NY)),VME_TOT)
!      call MPI_INST_PROF(dble(CR3(0,0,1:NY)),WME_TOT)
!      call MPI_INST_PROF(dble(CRTH(0,0,1:NY,1)),TH_TOT)
!      call MPI_INST_PROF(dble(CRTH(0,0,1:NY,2)),SAL_TOT) 
!      call MPI_INST_PROF(urms(1:NY),urms_tot)
!      call MPI_INST_PROF(vrms(1:NY),vrms_tot)
!      call MPI_INST_PROF(wrms(1:NY),wrms_tot)
!      call MPI_INST_PROF(shear(1:NY),shear_t)
!      call MPI_INST_PROF(dthdy(1:NY,1),dthdy_t)
!      call MPI_INST_PROF(dthdy(1:NY,2),dsdy_t)
!      call MPI_INST_PROF(thrms(1:NY,1),thrms_t)
!      call MPI_INST_PROF(thrms(1:NY,2),srms_t)
!      call MPI_INST_PROF(uv(1:NY),uv_t)
!      call MPI_INST_PROF(wv(1:NY),wv_t)
!      call MPI_INST_PROF(uw(1:NY),uw_t)
!      call MPI_INST_PROF(tke_mean(1:NY),tke_t)
!      call MPI_INST_PROF(tke_7(1:NY),diss_t)
!      call MPI_INST_PROF(tke_3(1:NY),prod_t)
!      call MPI_INST_PROF(bu_pro(1:NY),bpro_t)    
      if (MOVIE) then
! This file will contain a time history over a plane
         call SAVE_PLANE
!        call SAVE_PLANE_3D
      end if

 
      if (rank ==0)then
       tstep=TIME_STEP
        do j=1,NY*NPROCS
           write(105,*)urms_tot(j) 
          enddo
!       call NETCDF_OPEN_STATS_CHAN
       call NETCDF_WRITE_STATS_CHAN
!          do j=1,NY*NPROCS
!           write(105,*) ,tke_t(j)
!          enddo
      if (FINAL) then
! We are done with the simulation
! Close the NetCDF file
       CALL NETCDF_CLOSE_STATS_CHAN
      endif  

      open(300,file='output_folder/time_bulk_meltrate.txt',
     +    form='formatted',status='unknown',
     +    position='append')


      write(300,130) time, dble(CR2(0,0,1)),dble(CRTH(0,0,1,1)),
     &          dble(CRTH(0,0,1,2)),dthdz_mean(1), dthdz_mean(2),
     &          (C_sp_heat/L_heat)*(NU/Pr(1))*(1.0d0/0.92)
     &          *dthdz_mean(1),WME_TOT(NY*NPROCS),dwdy(2)

130   format (10f17.9)
      close(300)


      open(303,file='output_folder/instant_prof_tot.txt',
     +    form='formatted',status='unknown')
        do j=1,NY*NPROCS
         if(WME_TOT(j)/WME_TOT(NY*NPROCS) .ge. 0.90d0) then
!          vbl(j)=WME_TOT(j)/WME_TOT(NY*NPROCS)
!          if( vbl(j).ge.0.90d0)then
          vblf=GYF_TOT(j)
          goto 11
         endif
        enddo
11      write(6,*) 'the vbl is', vblf
      do j=1,NY*NPROCS
       write(303,566)GYF_TOT(j),UME_TOT(j),VME_TOT(j),WME_TOT(j),
     +  TH_TOT(j),SAL_TOT(j),urms_tot(j), vrms_tot(j), wrms_tot(j),
     +  shear(j),dthdy_t(j),dsdy_t(j),thrms_t(j),srms_t(j),tke_t(j),
     +  uv_t(j),wv_t(j),uw_t(j)  
!     +                SAL_TOT(j)
        write(306,*) thut_t(j), thvs_t(j)                      
      enddo
      close(303)
      endif

565    format(f12.5,4f19.8)
566    format(f12.5,18f19.8)



!      if (MOVIE) then
! This file will contain a time history over a plane
!         call SAVE_PLANE
!        call SAVE_PLANE_3D
!      end if

      do N=1,N_TH  
! Convert back to Fourier space
      call FFT_XZ_TO_FOURIER(TH(0,0,0,n),CTH(0,0,0,n),0,NY+1)

! End do over number of passive scalars, n
      end do
   
!       CALL tkebudget_chan_1
       


       if(rank.eq.0)then 
        open(44,file='wallnormal_pro.dat',status='unknown',
     +  form='formatted')!,position='append')
         write(44,*) TIME_STEP,TIME,DELTA_T
         do j=1,NY*NPROCS
            write(44,403) GYF_TOT(j),WME_TOT(j),UME_TOT(j),VME_TOT(j)
!     +  TH_TOT(j),SAL_TOT(j),urms_tot(j), vrms_tot(j), wrms_tot(j),
!     +  shear_t(j),dthdy_t(j),dsdy_t(j),thrms_t(j),srms_t(j),tke_t(j),
!     +  uv_t(j),wv_t(j),uw_t(j),diss_t(j),prod_t(j)
           write(663,407)!GYF_TOT(j),tke_t(j),prod_t(j),diss_t(j),
     +  bpro_t(j),thrms_t(j),dthdy_t(j),vrms_tot(j)!uv_t(j)
         enddo
!       call nc_write
       endif
403   format(4(F20.9,' '))
!       call nec_write
407    format(9(F20.9,''))
                          

!        IF (MOVIE) THEN
! This file will contain a single plane and is used in conjunction with
! the matlab script 'realtime_movie' to visualize data during
! simulation
!        open (76,file='temp.txt',status='unknown',form='formatted')
!        do J=1,NY
!          write(76,*) gyf(j)
!        end do
!        do I=0,NXM
!        do J=1,NY
!          write(76,*) U1(I,0,J)
!        end do
!        end do
!        close (76)
!        CALL SYSTEM('mv temp.txt ../post_process/matlab/latest_slice.txt
!     &')
!        END IF

! Write out the mean statistics at each time
      open(41,file=FNAME_TH,form='formatted',status='unknown',
     + position='append')
      write(41,*) TIME_STEP,TIME,DELTA_T,UBULK
      do n=1,N_TH 
      do j=1,NY
        write(41,402) j,GYF(J),dble(CTH(0,0,j,n))
     +      ,dthdy(j,n),thrms(j,n),thv(j,n),pe_diss(j,n)
      end do
      end do
      close(41)
402   format(I3,' ',6(F20.9,' '))

!      write(*,*) 'VERBOSITY: ',VERBOSITY
      if (VERBOSITY.gt.4) then 
      write(*,*) 'Outputting info for gnuplot...'
      open (unit=10, file="solution")
      do i=2,NXM
        do j=2,NYM
          write (10,*) i, j, U1(i,0,j)
        end do
        write (10,*) ""
      end do
      close (10)
      call system ('gnuplot <gnuplot.in') 
      end if


     


C Convert velocity back to Fourier space
      call fft_xz_to_fourier(U1(0,0,0),CU1(0,0,0),0,NY+1)
      call fft_xz_to_fourier(U2(0,0,0),CU2(0,0,0),0,NY+1)
      call fft_xz_to_fourier(U3(0,0,0),CU3(0,0,0),0,NY+1)
      call fft_xz_to_fourier(P(0,0,0),CP(0,0,0),0,NY+1)
  
      if (rank .eq. 0) then 
      write(*,*) 'done save_stats chan' 
      endif
      RETURN
      END

      subroutine tkebudget_chan
! NOte, it is important to only run this routine after complete R-K
!  time advancement since F1 is overwritten which is needed between R-K steps

      include 'header'

      integer i,j,k

! Compute the turbulent dissipation rate, epsilon=nu*<du_i/dx_j du_i/dx_j>
      do j=2,NYM
        epsilon(j)=0.
      end do
! Store du/dx in CS1
      do j=2,NYM
      do k=0,TNKZ
      do i=0,NKX
        CS1(i,k,j)=CIKX(i)*CR1(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=2,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Store dv/dx in CS1
      do j=2,NYM
      do k=0,TNKZ
      do i=0,NKX
        CS1(i,k,j)=CIKX(i)*(CR2(i,k,j)+CR2(i,k,j+1))/2.0
      end do
      end do
      end do
! Convert to physical space
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=2,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+0.5*(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Compute du/dy at GYF gridpoints, note remove mean
      do j=2,NYM
      do k=0,NZM
      do i=0,NXM
        F1(i,k,j)=((U1(i,k,j+1)-CR1(0,0,j+1))
     &      -(U1(i,k,j-1)-CR1(0,0,j-1)))/(GY(j)+GY(j+1))
      end do
      end do
      end do
      do j=2,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+0.5*(F1(i,k,j)**2.0)
! Cross term dvdx*dudy
        epsilon(j)=epsilon(j)+(S1(i,k,j)*F1(i,k,j))
      end do
      end do
      end do
! Store dw/dx in CS1
      do j=2,NYM
      do k=0,TNKZ
      do i=0,NKX
        CS1(i,k,j)=CIKX(i)*CR3(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=2,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+0.5*(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Compute du/dz at GYF gridpoints, note remove mean
! Store du/dz in CS1
      do j=2,NYM
      do k=0,TNKZ
      do i=0,NKX
        CF1(i,k,j)=CIKZ(k)*CR1(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xz_to_physical(CF1,F1,0,NY+1)
      do j=2,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+0.5*(F1(i,k,j)**2.0)
! Cross term dudz*dwdx
        epsilon(j)=epsilon(j)+S1(i,k,j)*F1(i,k,j)
      end do
      end do
      end do
! Compute dv/dy at GYF gridpoints, note remove mean
      do j=2,NYM
      do k=0,NZM
      do i=0,NXM
        S1(i,k,j)=((U2(i,k,j+1)-CR2(0,0,j+1))-(U2(i,k,j)-CR2(0,0,j)))
     &            /GYF(j)
      end do
      end do
      end do
      do j=2,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Compute dw/dy at GYF gridpoints, note remove mean
      do j=2,NYM
      do k=0,NZM
      do i=0,NXM
        S1(i,k,j)=((U3(i,k,j+1)-CR3(0,0,j+1))
     &      -(U3(i,k,j-1)-CR3(0,0,j-1)))/(GY(j)+GY(j+1))
      end do
      end do
      end do
      do j=2,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+0.5*(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Store dv/dz in CF1
      do j=2,NYM
      do k=0,TNKZ
      do i=0,NKX
        CF1(i,k,j)=CIKZ(k)*(CR2(i,k,j)+CR2(i,k,j+1))/2.0
      end do
      end do
      end do
! Convert to physical space
      call fft_xz_to_physical(CF1,F1,0,NY+1)
      do j=2,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+0.5*(F1(i,k,j)**2.0)
! Cross term dvdz*dwdy
        epsilon(j)=epsilon(j)+S1(i,k,j)*F1(i,k,j)
      end do
      end do
      end do
! Store dw/dz in CS1
      do j=2,NYM
      do k=0,TNKZ
      do i=0,NKX
        CS1(i,k,j)=CIKZ(k)*CR3(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=2,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+(S1(i,k,j)**2.0)
      end do
      end do
      end do
      do j=2,NYM
        epsilon(j)=NU*epsilon(j)/float(NX*NZ)
      end do


! Write out the bulk rms velocity
      !write(*,*) '<U_rms>: ',urms_b


! Write out the mean statistics at each time
      open(79,file='tke_epsi.txt',form='formatted',status='unknown')
c      write(45,*) TIME_STEP,TIME,DELTA_T
      do j=2,NYM
        write(79,401) j,GYF(J),epsilon(j)
      end do
401   format(I3,' ',2(F20.9,' '))


      return 
      end
       subroutine sponge

! This subroutine applies a sponge relaxation (Rayleigh damping) towards a
! specified background state for the temperature field
! The intention is to allow an open boundary
      include 'header'

! The following variables will store the background state

      real*8 b,co_b

      integer j

      NYC  = 45
      co_b = 41.d0
      b    = log(co_b)/(REAL(GYF(NY+1)-GYF(NYC)))

      DO J = 0, NY+1
        IF ( J .LT. NYC ) THEN
          SPONGE_SIGMA(J) = 0.d0
        ELSE
          SPONGE_SIGMA(J) = exp(b*(REAL(GYF(J)-GYF(NYC)))) -1.0
        ENDIF
      ENDDO

      return
      end       



      subroutine sponge_modi 

! This subroutine applies a sponge relaxation (Rayleigh damping) towards a
! specified background state for the temperature field
! The intention is to allow an open boundary
      include 'header'

! The following variables will store the background state

      real*8 b,co_b
      PARAMETER (co_b=1.4d0) 
      integer j

      b    = log(co_b)/(REAL(y_sp_max-y_sp_min))

      DO J = 0, NY+1
        IF ( J .LT. NYC ) THEN
          SPONGE_SIGMA(J) = 0.d0
        ELSE
          SPONGE_SIGMA(J) = exp(b*(REAL(GYF(J)-y_sp_min))) -1.0
        ENDIF
      ENDDO

      return
      end

      subroutine sponge_forcing

! This subroutine applies a sponge relaxation (Rayleigh damping) towards a
! specified background state for the temperature field
! The intention is to allow an open boundary
      include 'header' 
      INCLUDE "mpif.h"
      include 'header_mpi'

! The following variables will store the background state

      real*8 y_sp_forcing,y_sp_st,shift_sp,scal_sp,b,amp_f,gau_b
      PARAMETER (b=0.5d0, amp_f=0.45d0, gau_b=0.05d0) 
      integer i,j
      character*35 F_NUM

      I = RANK+1

      F_NUM ='grids/sponge_' 
     &        //CHAR(MOD(I,100)/10+48)
     &        //CHAR(MOD(I,10)+48) // '.txt'
      OPEN(111, file=F_NUM,form='formatted',status='unknown') 

      y_sp_forcing = 100000000.0 ;
      y_sp_st      = 100000000.0 ;           

      DO j=0,NY+1
       SPONGE_FORCE(j)=0.d0
       if ( (gyf(j) .gt. 100.d0 ) .and. (gyf(j) .lt. 102.d0) ) then
        y_sp_forcing = gyf(j) ;
       endif
      ENDDO

      CALL MPI_SECLT(y_sp_forcing)
      

      DO J = 0, NY+1
        SPONGE_FORCE(J) = -TANH(b*(gyf(j)-y_sp_forcing))
      ENDDO

      scal_sp  = SPONGE_FORCE(1) ;
      shift_sp = SPONGE_FORCE(NY+1) ;

      CALL  MPI_BCAST_REAL(scal_sp,1,1)
      CALL  MPI_BCAST_REAL_L(shift_sp,1,1)
      scal_sp = scal_sp-shift_sp;
       if (rank .eq.0)then
       write(6,*) 'Sponge_f', y_sp_forcing, rank,scal_sp,shift_sp 
       endif     
      DO J = 0, NY+1
        SPONGE_FORCE(J) = amp_f*(SPONGE_FORCE(J)-shift_sp)/scal_sp ;
!        WRITE(111,*) gyf(J),SPONGE_FORCE(J),SPONGE_SIGMA(j)
      ENDDO

      DO j=0,NY+1
       if ( (gyf(j) .gt. 16.d0 ) .and. (gyf(j) .lt. 17.d0) ) then
        y_sp_st = gyf(j) ;
       endif
      ENDDO

      CALL MPI_SECLT(y_sp_st)

      DO J = 0, NY+1
       if ( gyf(j) .le. y_sp_st ) then
        SPONGE_FORCE(J) = amp_f*exp(-gau_b*(gyf(J)-y_sp_st)**2. )
       else
        SPONGE_FORCE(J) = SPONGE_FORCE(J)                   ;
       endif
        WRITE(111,*) gyf(J),SPONGE_FORCE(J),SPONGE_SIGMA(j)
      ENDDO

      close(111)


      return
      end



      subroutine sponge_th(N)
! This subroutine applies a sponge relaxation (Rayleigh damping) towards a
! specified background state for the temperature field
! The intention is to allow an open boundary
      include 'header'

! The following variables will store the background state
      complex*8 TH_0(-1:NY+1)

      real*8 TH_TOP,b,co_b

      integer i,j,k,N



       IF (F_TYPE .EQ. 0) THEN
C     
        do k=0,TNKZ
         do i=0,NKX
          do j=0,NY+1
               CFTH(i,k,j,n)=CFTH(i,k,j,n)
     &                 -SPONGE_SIGMA(j)*(CTH(i,k,j,n)-0.d0)
          end do
         end do
        end do

      ELSEIF (F_TYPE .GT. 0) THEN
! Get mean at top of computational domain
!      TH_TOP=CTH(0,0,NYC,N)
! Do not let the value to float
       TH_TOP=GYF(NYC)

! Set the background state
      do j=NYC,NY+1
! Neumann
!        if ((TH_BC_UPPER(n).eq.1).or.(TH_BC_UPPER(n).eq.3)) then
          TH_0(j)=TH_TOP+(GYF(j)-GYF(NYC))*TH_BC_YMAX_C1(N)
!           TH_0(j)=GYF(j)
!        else if (TH_BC_UPPER(n).eq.0) then
! Dirichlet
!          TH_0(j)=TH_BC_UPPER_C1(n)
!        end if
c         write(199,*)j, TH_0(j), gyf(j), gyf(NYC)
      end do
! Add damping to R-K terms
c       do k=1,TNKZ
c         do i=1,NKX
c           do j=0,NY+1
c              CFTH(i,k,j,n)=CFTH(i,k,j,n)
c     &                 -SPONGE_SIGMA(j)*(CTH(i,k,j,n)-0.)
c           end do
c         end do
c      end do

! Damp the mean terms towards the background state
c      do j=0,NY+1
c        CFTH(0,0,j,n)=CFTH(0,0,j,n)
c     &        -SPONGE_SIGMA(j)*(CTH(0,0,j,n)-TH_0(j))
c        CFTH(1,0,j,n)=CFTH(1,0,j,n)
c     &        -SPONGE_SIGMA(j)*(CTH(1,0,j,n)-0.d0)
c        CFTH(0,1,j,n)=CFTH(0,1,j,n)
c     &        -SPONGE_SIGMA(j)*(CTH(0,1,j,n)-0.d0)
c      end do

      do k=0,TNKZ
         do i=0,NKX
           do j=0,NY+1
            If ((k.eq.0).and.(i.eq.0)) then
              CFTH(0,0,j,n)=CFTH(0,0,j,n)
     &        -SPONGE_SIGMA(j)*(CTH(0,0,j,n)-TH_0(j))
            else 
              CFTH(i,k,j,n)=CFTH(i,k,j,n)
     &                 -SPONGE_SIGMA(j)*(CTH(i,k,j,n)-0.)           
            endif
           end do
         end do
      end do

      endif

      return
      end


      subroutine sponge_vel
! This subroutine applies a sponge relaxation (Rayleigh damping) towards a
! specified background state
! The intention is to allow an open boundary
      include 'header'

! The following variables will store the background state
      real*8 U1_0(-1:NY+1), U2_0(0:NY+1) !, U3_0(-1:NY+1)

      integer i,j,k,N


! Set the background state
! Here, set the background to be geostrophic, with a linear temperature profile
      do j=0,NY+1
!        U1_0(j)=G_TAU/1.0
!        U3_0(j)=0.
      end do
      do j=0,NY+1
        U2_0(j)=0.
      end do
! Add damping function to explicit R-K
       do k=0,TNKZ
         do i=0,NKX
           if ((i.ne.0).or.(k.ne.0)) then
           do j=jstart,jend
             CF1(I,K,J)=CF1(I,K,J)-SPONGE_SIGMA(j)*(CU1(i,k,j)-0.d0)
             CF3(I,K,J)=CF3(I,K,J)-SPONGE_SIGMA(j)*(CU3(i,k,j)-0.d0)
           end do
           do j=1,NY
             CF2(I,K,J)=CF2(I,K,J)-
     &        0.5*(SPONGE_SIGMA(j)+SPONGE_SIGMA(j+1))*(CU2(i,k,j)-0.d0)
           end do
           end if
         end do
      end do
! Damp mean flow
      do j=jstart,jend
        CF1(0,0,j)=CF1(0,0,j)-SPONGE_SIGMA(j)*(CU1(0,0,j)-0.d0)
        CF3(0,0,j)=CF3(0,0,j)-SPONGE_SIGMA(j)*(CU3(0,0,j)-UBULK0)
      end do
      do j=1,NY
        CF2(0,0,j)=CF2(0,0,j)-SPONGE_SIGMA(j)*(CU2(0,0,j)-0.0d0)
      end do

      return
      end

      subroutine forcing_th
! This subroutine applies a sponge relaxation (Rayleigh damping) towards a
! specified background state
! The intention is to allow an open boundary
      include 'header'
      INCLUDE "mpif.h"
      include 'header_mpi'

! The following variables will store the background state
      real*8 b_ga,c_ga,f_force

      integer i,j,k,N
       c_ga = 150.0 ;
 
       IF (TIME .eq. 0.d0) THEN
        c_ga = 150.0 ;
        b_ga = 0.0 ;
        int_fc = 1 ;
        f_force = 0.0 ;
        TIME_fc_den = 0.d0 !3.0d0*c_ga ;
       ENDIF

       if (TIME .gt. TIME_fc_den) then
        f_force = 0.0 ;
        TIME_fc_den  = TIME_fc_den + 2.d0*pi/OMEGA0 !+ 1.0d0*c_ga ;
        int_fc  = int_fc + 1 ;
       else
        b_ga= 2*pi*real((int_fc - 1))/OMEGA0 ;
        f_force = exp(-(TIME-b_ga)**2.0d0/c_ga**2.0d0);
       endif

       if( (RK_STEP .eq. 3).and. (rank .eq. 1))THEN
       open(333,file='th_forcing.dat',form='formatted',status='unknown',
     + position='append')
       write(333,331) TIME,TIME/(2.0*pi/OMEGA0),f_force,b_ga,
     & TIME_fc_den, c_ga
       close(333)
331    format (2f16.6, f15.12, 3f16.6) 
       ENDIF

       do N=1,N_TH
       do k=0,TNKZ
         do i=0,NKX
          do j=0,NY+1
               CFTH(i,k,j,n)=CFTH(i,k,j,n)
     &             -f_force*SPONGE_FORCE(j)*(CTH(i,k,j,n)-0.d0)
          end do
         end do
        end do
       enddo  


       return
       end


      subroutine forcing_vel
! This subroutine applies a sponge relaxation (Rayleigh damping) towards a
! specified background state
! The intention is to allow an open boundary
      include 'header'
      INCLUDE "mpif.h"
      include 'header_mpi'

! The following variables will store the background state
      real*8 U1_0(-1:NY+1), U2_0(0:NY+1) ! U3_0(-1:NY+1)

      integer i,j,k,N

      real*8 A,L,s_1,s_2,s1z,s1zz,s2z,s2zz,
     &       t1,t1t
      real*8 p1,p2,p3,p4,p5,p6,p7,p8,p9,p10 
      
                  

        IF (CURVE_FITT) THEN
         t1  = cos(OMEGA0*time);
         IF (TIME .LE. 2.0) THEN
          A = 1.d0*TIME/(2.d0)
!         if (rank .eq. 0) write(6,*) 'Value of A', A
         ELSE
          A = 1.d0
         ENDIF

         do j=jstart,jend
           CF1(0,0,j)=CF1(0,0,j)
     &      - A*SPONGE_FORCE(J)*(CU1(0,0,j)-U3_0(j)*t1 ) 
         end do 

        ELSE
         do j=0,NY+1
          s_1 = exp(fact3*gyf(j)) ;
          s_2 = sin(fact4*gyf(j)) ;
          t1  = cos(OMEGA0*time);
          U3_0(j) = (u_0/0.47)*s_1*s_2*t1 ;

          s1z =fact3*exp(fact3*gyf(j)) ;
          s1zz=fact3**2.0*exp(fact3*gyf(j)) ;
          s2z=fact4*cos(fact4*gyf(j)) ;
          s2zz=-fact4**2.0*sin(fact4*gyf(j)) ;
          t1t=-omega0*sin(omega0*time);

          U1_0(j)=NU*u_0/0.47*t1*(s1zz*s_2+2.0*s1z*s2z+s_1*s2zz)
!     &         - u_0/0.47*s_1*s_2*t1t 
!     &        - RI_FINAL*CTH(0,0,j,1)*sin( ANG_BETA)  ;

         end do

         IF (TIME .LE. 2.0) THEN
          A = 1.d0*TIME/(2.d0)
!          if (rank .eq. 0) write(6,*) 'Value of A', A
         ELSE
           A = 1.d0
         ENDIF

         do j=jstart,jend
          CF1(0,0,j)=CF1(0,0,j) 
     &     - A*SPONGE_FORCE(J)*(CU1(0,0,j)-U3_0(j)) -U1_0(j)
         end do
      ENDIF

      return
      end

       subroutine tkebudget_chan_1
! NOte, it is important to only run this routine after complete R-K
!  time advancement since F1 is overwritten which is needed between R-K steps
! This subroutine should be called in SAVE_STATS_CHAN after computing
! plane averaged statistics, with the velocity in physical space, and 
! CRi containing the velocity in Fourier space

      include 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'
 
      integer i,j,k,n
      CHARACTER*35 FNAME
! Define working arrays
      real*8 tke_1(0:NY+1)
      real*8 tke_2(0:NY+1)
      real*8 tke_2_1(0:NY+1)
      real*8 tke_2_2(0:NY+1)
      real*8 tke_2_3(0:NY+1)
!      real*8 tke_3(0:NY+1)
      real*8 tke_3_1(0:NY+1)
      real*8 tke_3_2(0:NY+1)
      real*8 tke_3_3(0:NY+1)
      real*8 tke_4(0:NY+1)
      real*8 tke_5(0:NY+1,1:N_TH)
      real*8 tke_6_1(0:NY+1)
      real*8 tke_6_1_1(0:NY+1)
      real*8 tke_6_1_2(0:NY+1)
      real*8 tke_6_1_3(0:NY+1)
      real*8 tke_6_2(0:NY+1)
      real*8 tke_6_2_1(0:NY+1)
      real*8 tke_6_2_2(0:NY+1)
      real*8 tke_6_2_3(0:NY+1)
!      real*8 tke_7(0:NY+1)
      real*8 tke_7_1(0:NY+1)
      real*8 tke_7_2(0:NY+1)
      real*8 tke_7_3(0:NY+1)
      real*8 tke_8_1(0:NY+1)
      real*8 tke_8_2(0:NY+1)
      real*8 tke_8_3(0:NY+1)
      real*8 S1_mean(0:NY+1)
      real*8 p_mean(0:NY+1)
      real*8 F2_mean(0:NY+1)
      real*8 transport(0:NY+1)


      I = RANK+1

           FNAME ='results/tke_'
     &        //CHAR(MOD(I,100)/10+48)
     &        //CHAR(MOD(I,10)+48) // '.txt'


! Preparation work:
! Store the U2 velocity interpolated to the GYF grid in F2
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            F2(i,k,j)=0.5d0*(U2(i,k,j+1)+U2(i,k,j))
          end do
        end do
      end do       
      do j=1,NY
        F2_mean(j)=0.5d0*(dble(CR2(0,0,J+1))+dble(CR2(0,0,J)))
      end do

      do j=1,NY+1
! tke_mean defined at GY points
        tke_mean_old(j)=tke_mean(j)
c        tke_mean(j)=0.5d0*((0.5*(urms(j-1)+urms(j)))**2.d0
c     :       +vrms(j)**2.d0
c     :       +(0.5*(wrms(j-1)+wrms(j)))**2.d0) 

         tke_mean(j) = 0.5*((DYF(J-1)*urms(j)**2.d0
     &                 + DYF(j)*urms(j-1)**2.d0)/(2.d0*DY(J))
     &                 + vrms(j)**2.d0
     &                 + (DYF(J-1)*wrms(j)**2.d0
     &                 + DYF(j)*wrms(j-1)**2.d0)/(2.d0*DY(J)))

        IF (TIME_STEP > 1 ) THEN
        tke_1(j)=(tke_mean(j)-tke_mean_old(j))/(TIME-TIME_old)
        ELSE
        tke_1(j)= 0.d0
        ENDIF 
      end do
      tke2d(:,:)=0.0d0
      do j=1,NY+1
       do i=1,NX
! tke_mean2d defined at GY points
         tke2d(i,j) = 0.5/(NZ)*((DYF(J-1)*urms2d(i,j)**2.d0
     &                 + DYF(j)*urms2d(i,j-1)**2.d0)/(2.d0*DY(J))
     &                 + vrms2d(i,j)**2.d0
     &                 + (DYF(J-1)*wrms2d(i,j)**2.d0
     &                 + DYF(j)*wrms2d(i,j-1)**2.d0)/(2.d0*DY(J)))

        end do
       enddo
       do j=1,NY
         write(414,*) tke2d(100,j),DYF(j), wrms2d(100,j)
       enddo
        
       time_old=TIME

      do j=2,NY
!        tke_2(j)=-dble(CR2(0,0,J))*(tke_mean(j+1)-tke_mean(j-1))
!     &        /(GYF(J+1)-GYF(J-1))
        tke_2(j)=-dble(CR2(0,0,J))*(tke_mean(j+1)-tke_mean(j-1))
     &        /(GY(J+1)-GY(j-1))

      end do


! Get the individual terms
! Use S1 as a working variable
! U1*U1 
      do j=1,NY
        do k=0,TNKZ
          do i=0,NKX
            CS1(i,k,j)=CIKX(I)*CR1(i,k,j)
          end do
        end do
      end do
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=(U1(i,k,j)-dble(CR1(0,0,j)))
     &                      *S1(i,k,j)*dble(CR1(0,0,j))
          end do
        end do
      end do
      do j=1,NY
        tke_2_1(j)=-1.d0*SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do       
! U2*U1
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=(U1(i,k,j)-dble(CR1(0,0,j))
     &             *(0.5d0*(U1(i,k,j+1)-dble(CR1(0,0,j+1))
     &                     +U1(i,k,j)-dble(CR1(0,0,j)))
     &               -0.5d0*(U1(i,k,j)-dble(CR1(0,0,j))
     &                      +U1(i,k,j-1)-dble(CR1(0,0,j-1))))/DYF(j))
     &             *0.5d0*(CR2(0,0,j+1)+CR2(0,0,j))
          end do
        end do
      end do
      do j=1,NY
        tke_2_1(j)=tke_2_1(j)-SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do       
! U3*U1 
      do j=1,NY
        do k=0,TNKZ
          do i=0,NKX
            CS1(i,k,j)=CIKZ(K)*CR1(i,k,j)
          end do
        end do
      end do
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=(U1(i,k,j)-CR1(0,0,j))*S1(i,k,j)*CR3(0,0,j)
          end do
        end do
      end do
      do j=1,NY
        tke_2_1(j)=tke_2_1(j)-SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do
! U1*U2
      do j=2,NY       
        do k=0,TNKZ
          do i=0,NKX
            CS1(i,k,j)=CIKX(I)*CR2(i,k,j)
          end do
        end do
      end do

      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=2,NY 
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=(U2(i,k,j)-dble(CR2(0,0,j)))*S1(i,k,j)
     &            *(DYF(J)*dble(CR1(0,0,j))+DYF(J-1)*dble(CR1(0,0,j-1)))
     &             /(2.d0*DY(J))
          end do
        end do
      end do
      do j=2,NY
        tke_2_2(j)=-1.d0*SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do
! U2*U2
      do j=2,NY
        do k=0,NZM  
          do i=0,NXM
            S1(i,k,j)=dble(CR2(0,0,j))*(U2(i,k,j)-dble(CR2(0,0,j)))
     &              *(0.5d0*(U2(i,k,j+1)+U2(i,k,j))
     &               -0.5d0*(U2(i,k,j)+U2(i,k,j-1)))/DY(J)
          end do
        end do
      end do
      do j=2,NY
        tke_2_2(j)=tke_2_2(j)-SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do
! U2*U3
      do j=2,NY
        do k=0,TNKZ
          do i=0,NKX
            CS1(i,k,j)=CIKZ(K)*CR2(i,k,j)
          end do
        end do
      end do
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=2,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=(U2(i,k,j)-dble(CR2(0,0,j)))*S1(i,k,j)
     &            *(DYF(J)*dble(CR3(0,0,j))+DYF(J-1)*dble(CR3(0,0,j-1)))
     &             /(2.d0*DY(J))
          end do
        end do
      end do
      do j=2,NY
        tke_2_2(j)=tke_2_2(j)-SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do
! U3*U1 
      do j=1,NY
        do k=0,TNKZ
          do i=0,NKX
            CS1(i,k,j)=CIKX(I)*CR3(i,k,j)
          end do
        end do
      end do
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=(U3(i,k,j)-dble(CR3(0,0,j)))
     &         *S1(i,k,j)*dble(CR3(0,0,j))
          end do
        end do
      end do
      do j=1,NY
        tke_2_3(j)=-1.d0*SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do
! U2*U3
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=(U3(i,k,j)-dble(CR3(0,0,j))
     &             *(0.5d0*(U3(i,k,j+1)-dble(CR3(0,0,j+1))
     &                     +U3(i,k,j)-dble(CR3(0,0,j)))
     &               -0.5d0*(U3(i,k,j)-dble(CR3(0,0,j))
     &                      +U3(i,k,j-1)-dble(CR3(0,0,j-1))))/DYF(j))
     &             *0.5d0*(CR2(0,0,j+1)+CR2(0,0,j))
          end do
        end do
      end do
      do j=1,NY
        tke_2_3(j)=tke_2_3(j)-SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do
! U3*U3 
      do j=1,NY
        do k=0,TNKZ
          do i=0,NKX
            CS1(i,k,j)=CIKZ(K)*CR3(i,k,j)
          end do
        end do
      end do
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=(U3(i,k,j)-CR3(0,0,j))*S1(i,k,j)*CR3(0,0,j)
          end do
        end do
      end do
      do j=1,NY
        tke_2_3(j)=tke_2_3(j)-SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do


! Get the production at GY points
      do j=2,NY
!        tke_3(j)=-uv(j)*dudy(j)-wv(j)*dwdy(j)
         tke_3(j)=-uv(j)*(dble(CR1(0,0,j))-dble(CR1(0,0,j-1)))/DY(j)
     &           -wv(j)*(dble(CR3(0,0,j))-dble(CR3(0,0,j-1)))/DY(j)
      end do
        tke_3(NY+1)=tke_3(NY)
        tke_3(1)=tke_3(2)

! Get the components of the production
! Use S1 as a working variable
! U1*U1
!      do j=1,NY
!       do k=0,TNKZ
!          do i=0,NKX
!            CS1(i,k,j)=CIKX(I)*CR1(i,k,j)
!          end do
!        end do
!      end do
!      call fft_xz_to_physical(CS1,S1,0,NY+1)
!      do j=1,NY
!        do k=0,NZM
!          do i=0,NXM
!            S1(i,k,j)=(U1(i,k,j)-dble(CR1(0,0,j)))
!     &                      *S1(i,k,j)*dble(CR1(0,0,j))
!          end do
!        end do
!      end do
!      do j=2,NY
!        tke_3_1(j)=-1.d0*SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
!      end do

      tke_3_1(1)=0.d0
! U2*U1
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
c            S1(i,k,j)=
c     &     (    0.5d0*(dble(CR1(0,0,j+1))+dble(CR1(0,0,j)))
c     &       - 0.5d0*(dble(CR1(0,0,j))+dble(CR1(0,0,j-1)))  ) /(DYF(j))
c            S1(i,k,j)=S1(i,k,j)*(U1(i,k,j)-dble(CR1(0,0,j)))
c     &                *(U2(i,k,j+1)-dble(CR2(0,0,j+1)))
             S1(i,k,j)=
     &              ( (DYF(J)*dble(CR1(0,0,j+1))
     &             + DYF(j+1)*dble(CR1(0,0,j)))/(2.d0*DY(J+1))
     &             - (DYF(J-1)*dble(CR1(0,0,j))
     &             + DYF(J)*dble(CR1(0,0,j-1)))/(2.d0*DY(J)) ) /(DYF(j))

            S1(i,k,j)=S1(i,k,j)*(U1(i,k,j)-dble(CR1(0,0,j)))
     &                *(0.5d0*(U2(i,k,j+1)+U2(i,k,j))
     &                - 0.5d0*(dble(CR2(0,0,j+1))+dble(CR2(0,0,j))))
 

          end do
        end do
      end do

      do j=1,NY
        tke_3_1(j)=-SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do
! U3*U1
!      do j=1,NY
!        do k=0,TNKZ
!          do i=0,NKX
!            CS1(i,k,j)=CIKZ(K)*CR3(i,k,j)
!          end do
!        end do
!      end do
!      call fft_xz_to_physical(CS1,S1,0,NY+1)
!      do j=1,NY
!        do k=0,NZM
!          do i=0,NXM
!            S1(i,k,j)=(U1(i,k,j)-CR1(0,0,j))*S1(i,k,j)*dble(CR1(0,0,j))
!          end do
!        end do
!      end do
!      do j=2,NY
!        tke_3_1(j)=tke_3_1(j)-SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
!      end do
! U1*U2
!      do j=2,NY
!        do k=0,TNKZ
!          do i=0,NKX
!            CS1(i,k,j)=CIKX(I)*(DYF(J)*CR1(i,k,j)+DYF(J-1)*CR1(i,k,j-1))
!     &             /(2.d0*DY(J))
!          end do
!        end do
!      end do
!      call fft_xz_to_physical(CS1,S1,0,NY+1)
!      do j=2,NY
!        do k=0,NZM
!          do i=0,NXM
!            S1(i,k,j)=(U2(i,k,j)-dble(CR2(0,0,j)))*S1(i,k,j)
!     &              *dble(CR2(0,0,j))
!          end do
!        end do
!      end do
!      do j=2,NY
!        tke_3_2(j)=-1.d0*SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
!      end do
! U2*U2
      do j=2,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=(U2(i,k,j)-dble(CR2(0,0,j)))**2.0
     &              *(0.5d0*(U2(i,k,j+1)+U2(i,k,j))
     &               -0.5d0*(U2(i,k,j)+U2(i,k,j-1)))/DY(J)
          end do
        end do
      end do
      do j=2,NY
        tke_3_2(j)=-SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do
! U3*U2
!      do j=2,NY
!        do k=0,TNKZ
!          do i=0,NKX
!            CS1(i,k,j)=CIKZ(K)*(DYF(J)*CR3(i,k,j)+DYF(J-1)*CR3(i,k,j-1))
!     &             /(2.d0*DY(J))
!          end do
!        end do
!      end do
!      call fft_xz_to_physical(CS1,S1,0,NY+1)
!      do j=2,NY
!        do k=0,NZM
!          do i=0,NXM
!            S1(i,k,j)=(U2(i,k,j)-dble(CR2(0,0,j)))*S1(i,k,j)
!     &              *dble(CR2(0,0,j))
!          end do
!        end do
!      end do
!      do j=2,NY
!        tke_3_2(j)=tke_3_2(j)-SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
!      end do
! U3*U1
!      do j=1,NY
!        do k=0,TNKZ
!          do i=0,NKX
!            CS1(i,k,j)=CIKX(I)*CR1(i,k,j)
!          end do
!        end do
!      end do
!      call fft_xz_to_physical(CS1,S1,0,NY+1)
!      do j=1,NY
!        do k=0,NZM
!          do i=0,NXM
!            S1(i,k,j)=(U3(i,k,j)-dble(CR3(0,0,j)))
!     &         *S1(i,k,j)*dble(CR3(0,0,j))
!          end do
!        end do
!      end do
!      do j=2,NY
!        tke_3_3(j)=-1.d0*SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
!      end do
      tke_3_3(1)=0.d0
! U2*U1
c      do j=1,NY
c        do k=0,NZM
c          do i=0,NXM
c            S1(i,k,j)=
c     &     ( 0.5d0*(dble(CR3(0,0,j+1))+dble(CR3(0,0,j)))
c     &       - 0.5d0*(dble(CR3(0,0,j))+dble(CR3(0,0,j-1)))  ) /(DYF(j))
c            S1(i,k,j)=S1(i,k,j)*(U2(i,k,j+1)-dble(CR2(0,0,j+1)))
c     &                *(U3(i,k,j)-dble(CR3(0,0,j)))
c          end do
c        end do
c      end do

       do j=1,NY
         do k=0,NZM
          do i=0,NXM 

             S1(i,k,j) = ( (DYF(J)*dble(CR3(0,0,j+1))
     &             + DYF(j+1)*dble(CR3(0,0,j)))/(2.d0*DY(J+1))
     &             - (DYF(J-1)*dble(CR3(0,0,j))
     &             + DYF(J)*dble(CR3(0,0,j-1)))/(2.d0*DY(J)) ) /(DYF(j))

             S1(i,k,j) = S1(i,k,j)*(U3(i,k,j)-dble(CR3(0,0,j)))
     &                   *(0.5d0*(U2(i,k,j+1)+U2(i,k,j))
     &                - 0.5d0*(dble(CR2(0,0,j+1))+dble(CR2(0,0,j))))

             end do
         end do
      end do


      do j=1,NY
        tke_3_3(j)=-SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do
! U3*U1
!      !do j=1,NY
!        do k=0,TNKZ
!          do i=0,NKX
!            CS1(i,k,j)=CIKZ(K)*CR3(i,k,j)
!          end do
!        end do
!      end do
!     call fft_xz_to_physical(CS1,S1,0,NY+1)
!      do j=1,NY
!        do k=0,NZM
!          do i=0,NXM
!            S1(i,k,j)=(U3(i,k,j)-CR3(0,0,j))*S1(i,k,j)*dble(CR3(0,0,j))
!          end do
!        end do
!      end do
!      do j=2,NY
!        tke_3_3(j)=tke_3_3(j)-SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
!      end do
      bp2d(:,:)=0.0d0
      do j=2,NY
        tke_4(j)=(NU/DY(J))*((tke_mean(j+1)-tke_mean(j))/DYF(J)
     &                        -(tke_mean(j)-tke_mean(j-1))/DYF(J-1))
      end do 

      do n=1,N_TH
        do j=2,NY
          tke_5(j,n)= RI_TAU(n)*(thu(j,n)*sin(ANG_BETA)-
     &                 thv(j,n)*cos(ANG_BETA))
        end do
        do i=1,NX
         do j=1,NY
          bp2d(i,j)=bp2d(i,j)+ RI_TAU(n)*(thu2d(i,j,n)*sin(ANG_BETA)-
     &                 thv2d(i,j,n)*cos(ANG_BETA))
         enddo
        enddo
      end do 
        do i=1,NX
           do j=1,NY
            bp2d(i,j)=bp2d(i,j)/float(NZ)
           enddo
            bp2d(i,NY+1)=bp2d(i,NY)
        enddo
      bu_pro(:)=0.0d0
      do j=2,NY
        do n=1,N_TH
          bu_pro(j)=bu_pro(j)+tke_5(j,n)
        enddo
      enddo
      bu_pro(1)=bu_pro(2)
      bu_pro(nY+1)=bu_pro(NY)  
! Construct the resolved turbulent transport terms 
! Convert the pressure to physical space for computation of the
! pressure transport
! Calcuate transport at GYF points since we will differentiate it d/dx2
! Already pressure field in fourier space is stored in CFTH
c      do j=2,NY
c        p_mean(j)=dble(CF1(0,0,j))
c      end do
      do j=1,NY
      transport(j)=0.d0
        do k=0,NZM
          do i=0,NZM
            transport(j)=transport(j)
! Pressure Transport term:
     &        +(F2(I,K,J)-F2_mean(j))*(P(I,K,J)-dble(CFTH(0,0,J,1)))
          end do
        end do
        transport(j)=transport(j)/dble(NX*NZ)
      end do
      do j=2,NY
        tke_6_1(j)=-(transport(j)-transport(j-1))/DY(J) 
      end do

! Get individual terms 
! 1

      do j=1,NY
        do k=0,TNKZ
          do i=0,NKX
            CS1(i,k,j)=CIKX(I)*CFTH(i,k,j,1)
          end do
        end do
      end do
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=S1(i,k,j)*(U1(i,k,j)-dble(CR1(0,0,j)))
          end do
        end do
      end do
      do j=1,NY
        tke_6_1_1(j)=-1.d0*SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do
! 3
      do j=1,NY
        do k=0,TNKZ
          do i=0,NKX
            CS1(i,k,j)=CIKZ(K)*CFTH(i,k,j,1)
          end do
        end do
      end do
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=S1(i,k,j)*(U3(i,k,j)-dble(CR3(0,0,j)))
          end do
        end do
      end do
      do j=1,NY
        tke_6_1_3(j)=-1.d0*SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do
! 2
c      do j=1,NY
c        p_mean(j)=dble(CP(0,0,J))
c      end do
c      call fft_xz_to_physical(CP,P,0,NY+1)

      do j=2,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=(U2(i,k,j)-dble(CR2(0,0,j)))
     &        *((P(i,k,j)-dble(CFTH(0,0,j,1)))-(P(i,k,j-1)
     &            -dble(CFTH(0,0,j-1,1))))/DY(J)
          end do
        end do
      end do 
      do j=2,NY
        tke_6_1_2(j)=-1.d0*SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ) 
      end do

! Turbulent transport terms:
      do j=1,NY
      transport(j)=0.d0
        do k=0,NZM
          do i=0,NZM
            transport(j)=transport(j)
! U1^2*U2
     &     +0.5d0*(U1(I,K,J)-dble(CR1(0,0,J)))**2.d0
     &           *(F2(I,K,J)-F2_mean(j))
! U2^3
     &        +0.5d0*(F2(I,K,J)-F2_mean(J))**3.d0
! U3^2*U2
     &     +0.5d0*(U3(I,K,J)-dble(CR3(0,0,J)))**2.d0
     &           *(F2(I,K,J)-F2_mean(j))
          end do
        end do
        transport(j)=transport(j)/dble(NX*NZ)
      end do
! We have transport defined at GYF points
! Now, the vertical derivative of the transport term:
      do j=2,NY
        tke_6_2(j)=-(transport(j)-transport(j-1))/DY(J) 
      end do

! Get the individual components of the turbulent transport term
! 1,1
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=(U1(i,k,j)-dble(CR1(0,0,j)))**2.d0
          end do
        end do
      end do
      call fft_xz_to_fourier(S1,CS1,0,NY+1)
      do j=1,NY
        do k=0,TNKZ
          do i=0,NKX
            CS1(i,k,j)=CIKX(I)*CS1(i,k,j)      
          end do
        end do
      end do
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=S1(i,k,j)*(U1(i,k,j)-dble(CR1(0,0,j)))
          end do
        end do 
      end do
      do j=1,NY
        tke_6_2_1(j)=-1.d0*SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do
! 1,2
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=
     &             (   0.5d0*((U1(i,k,j+1)-dble(CR1(0,0,j+1)))
     &                      +(U1(i,k,j)-dble(CR1(0,0,j))))
     &                     *(U2(i,k,j+1)-dble(CR2(0,0,j+1)))
     &               - 0.5d0*((U1(i,k,j)-dble(CR1(0,0,j)))
     &                      +(U1(i,k,j-1)-dble(CR1(0,0,j-1))))
     &                     *(U2(i,k,j)-dble(CR2(0,0,j)))  ) / DYF(j)
            S1(i,k,j)=S1(i,k,j)*(U1(i,k,j)-dble(CR1(0,0,j))) 
          end do
        end do
      end do 
      do j=1,NY
        tke_6_2_1(j)=tke_6_2_1(j)-1.d0*SUM(S1(0:NXM,0:NZM,j))
     &                                   /dble(NX*NZ)
      end do
! 1,3
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=(U1(i,k,j)-dble(CR1(0,0,j)))
     &               *(U3(i,k,j)-dble(CR3(0,0,j)))
          end do
        end do
      end do
      call fft_xz_to_fourier(S1,CS1,0,NY+1)
      do j=1,NY
        do k=0,TNKZ
          do i=0,NKX
            CS1(i,k,j)=CIKZ(K)*CS1(i,k,j)
          end do
        end do
      end do
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=S1(i,k,j)*(U1(i,k,j)-dble(CR1(0,0,j)))
          end do
        end do
      end do
      do j=1,NY
        tke_6_2_1(j)=tke_6_2_1(j)-1.d0*SUM(S1(0:NXM,0:NZM,j))
     &                                 /dble(NX*NZ)
      end do      
! 2,1
      do j=2,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=(U2(i,k,j)-dble(CR2(0,0,j)))
     &       *(DYF(J)*(U1(i,k,j)-dble(CR1(0,0,j)))
     &        +DYF(J-1)*(U1(i,k,j-1)-dble(CR1(0,0,j-1))))/(2.d0*DY(J))
          end do
        end do
      end do
      call fft_xz_to_fourier(S1,CS1,0,NY+1)
      do j=2,NY
        do k=0,TNKZ
          do i=0,NKX
            CS1(i,k,j)=CIKX(I)*CS1(i,k,j)
          end do
        end do
      end do
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=2,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=S1(i,k,j)*(U2(i,k,j)-dble(CR2(0,0,j)))
          end do
        end do
      end do
      do j=1,NY
        tke_6_2_2(j)=-1.d0*SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do    
! 2,2 
      do j=2,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=
     &             (   (0.5d0*((U2(i,k,j+1)-dble(CR2(0,0,j+1)))
     &                      +(U2(i,k,j)-dble(CR2(0,0,j)))))**2.d0
     &              -  (0.5d0*((U2(i,k,j)-dble(CR2(0,0,j)))
     &                      +(U2(i,k,j-1)-dble(CR2(0,0,j-1)))))**2.d0
     &                      ) / DYF(j)
            S1(i,k,j)=S1(i,k,j)*(U2(i,k,j)-dble(CR2(0,0,j)))
          end do
        end do
      end do
      do j=2,NY
        tke_6_2_2(j)=tke_6_2_2(j)-1.d0*SUM(S1(0:NXM,0:NZM,j))
     &                                   /dble(NX*NZ)
      end do
! 2,3 
      do j=2,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=(U2(i,k,j)-dble(CR2(0,0,j)))
     &       *(DYF(J)*(U3(i,k,j)-dble(CR3(0,0,j)))
     &        +DYF(J-1)*(U3(i,k,j-1)-dble(CR3(0,0,j-1))))/(2.d0*DY(J))
          end do
        end do
      end do
      call fft_xz_to_fourier(S1,CS1,0,NY+1)
      do j=2,NY
        do k=0,TNKZ
          do i=0,NKX
            CS1(i,k,j)=CIKZ(K)*CS1(i,k,j)
          end do
        end do
      end do
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=2,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=S1(i,k,j)*(U2(i,k,j)-dble(CR2(0,0,j)))
          end do
        end do
      end do
      do j=1,NY
        tke_6_2_2(j)=tke_6_2_2(j)-1.d0*SUM(S1(0:NXM,0:NZM,j))
     &                                   /dble(NX*NZ)
      end do
! 3,1
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=(U3(i,k,j)-dble(CR3(0,0,j)))
     &                *(U1(i,k,j)-dble(CR1(0,0,j)))
          end do
        end do
      end do
      call fft_xz_to_fourier(S1,CS1,0,NY+1)
      do j=1,NY
        do k=0,TNKZ
          do i=0,NKX
            CS1(i,k,j)=CIKX(I)*CS1(i,k,j)
          end do
        end do
      end do
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=S1(i,k,j)*(U3(i,k,j)-dble(CR3(0,0,j)))
          end do
        end do
      end do
      do j=1,NY
        tke_6_2_3(j)=-1.d0*SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do
! 3,2
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=
     &             (   0.5d0*((U3(i,k,j+1)-dble(CR3(0,0,j+1)))
     &                      +(U3(i,k,j)-dble(CR3(0,0,j))))
     &                     *(U2(i,k,j+1)-dble(CR2(0,0,j+1)))
     &               - 0.5d0*((U3(i,k,j)-dble(CR3(0,0,j)))
     &                      +(U3(i,k,j-1)-dble(CR3(0,0,j-1))))
     &                     *(U2(i,k,j)-dble(CR2(0,0,j)))  ) / DYF(j)
            S1(i,k,j)=S1(i,k,j)*(U3(i,k,j)-dble(CR3(0,0,j)))
          end do
        end do
      end do
      do j=1,NY
        tke_6_2_3(j)=tke_6_2_3(j)-1.d0*SUM(S1(0:NXM,0:NZM,j))
     &                                   /dble(NX*NZ)
      end do
! 3,3
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=(U3(i,k,j)-dble(CR3(0,0,j)))
     &               *(U3(i,k,j)-dble(CR3(0,0,j)))
          end do
        end do
      end do
      call fft_xz_to_fourier(S1,CS1,0,NY+1)
      do j=1,NY
        do k=0,TNKZ
          do i=0,NKX
            CS1(i,k,j)=CIKZ(K)*CS1(i,k,j)
          end do
        end do
      end do
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=S1(i,k,j)*(U3(i,k,j)-dble(CR3(0,0,j)))
          end do
        end do
      end do
      do j=1,NY
        tke_6_2_3(j)=tke_6_2_3(j)-1.d0*SUM(S1(0:NXM,0:NZM,j))
     &                                     /dble(NX*NZ)
      end do

 
! Compute the turbulent dissipation rate, epsilon=nu*<du_i/dx_j du_i/dx_j>
      do j=1,NY
        epsilon(j)=0.
      end do
! Store du/dx in CS1
      do j=1,NY
      do k=0,TNKZ
      do i=0,NKX
        CS1(i,k,j)=CIKX(i)*CR1(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=2,NY
      do k=0,NZM
      do i=0,NXM
!        epsilon(j)=epsilon(j)+(S1(i,k,j)**2.0)
         epsilon(j)=epsilon(j)+(DYF(j-1)*S1(i,k,j)**2.d0
     &             +DYF(j)*S1(i,k,j-1)**2.d0)/(2.d0*DY(j))
!         if(k==10) then
         ep2d(i,j)=ep2d(i,j)+(DYF(j-1)*S1(i,k,j)**2.d0
     &             +DYF(j)*S1(i,k,j-1)**2.d0)/(2.d0*DY(j))
!         endif
      end do
      end do
      end do
      write(66,*) ep2d(200,5)
! Store dv/dx in CS1
      do j=1,NY
      do k=0,TNKZ
      do i=0,NKX
!        CS1(i,k,j)=CIKX(i)*(CR2(i,k,j)+CR2(i,k,j+1))/2.0
        CS1(i,k,j)=CIKX(i)*CR2(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=2,NY
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+(S1(i,k,j)**2.0)
!        if(k==10) then
        ep2d(i,j)=ep2d(i,j)+(S1(i,k,j)**2.0) 
!        endif
      end do
      end do
      end do
      write(66,*) ep2d(200,5)
! Compute du/dy at GYF gridpoints, note remove mean
      do j=2,NY
      do k=0,NZM
      do i=0,NXM
!        F1(i,k,j)=((U1(i,k,j+1)-CR1(0,0,j+1))
!     &      -(U1(i,k,j-1)-CR1(0,0,j-1)))/(GYF(j+1)+GYF(j-1))
         F1(i,k,j)=((U1(i,k,j)-dble(CR1(0,0,j)))
     &          -(U1(i,k,j-1)-dble(CR1(0,0,j-1))))
     &             /DY(j)

      end do
      end do
      end do
      do j=2,NY
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+(F1(i,k,j)**2.0)
!        if(k==10) then
         ep2d(i,j)=ep2d(i,j)+(F1(i,k,j)**2.0)
!        endif 
      end do
      end do
      end do
      write(66,*) ep2d(200,5)
! Store dw/dx in CS1
      do j=1,NY
      do k=0,TNKZ
      do i=0,NKX
        CS1(i,k,j)=CIKX(i)*CR3(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=2,NY
      do k=0,NZM
      do i=0,NXM
!        epsilon(j)=epsilon(j)+(S1(i,k,j)**2.0)
         epsilon(j)=epsilon(j)+(DYF(j-1)*S1(i,k,j)**2.d0
     &             +DYF(j)*S1(i,k,j-1)**2.d0)/(2.d0*DY(j))
!         if(k==10) then
          ep2d(i,j)=ep2d(i,j)+(DYF(j-1)*S1(i,k,j)**2.d0
     &             +DYF(j)*S1(i,k,j-1)**2.d0)/(2.d0*DY(j))
!         endif
      end do
      end do
      end do
      write(66,*) ep2d(200,5)
! Compute du/dz at GYF gridpoints
! Store du/dz in CF1
      do j=1,NY
      do k=0,TNKZ
      do i=0,NKX
        CF1(i,k,j)=CIKZ(k)*CR1(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xz_to_physical(CF1,F1,0,NY+1)
      do j=2,NY
      do k=0,NZM
      do i=0,NXM
!        epsilon(j)=epsilon(j)+(F1(i,k,j)**2.0)
         epsilon(j)=epsilon(j)+(DYF(j-1)*F1(i,k,j)**2.d0
     &             +DYF(j)*F1(i,k,j-1)**2.d0)/(2.d0*DY(j))
!          if(k==10) then
         ep2d(i,j)=ep2d(i,j)+(DYF(j-1)*F1(i,k,j)**2.d0
     &             +DYF(j)*F1(i,k,j-1)**2.d0)/(2.d0*DY(j))
!          endif
      end do
      end do
      end do
       write(66,*) ep2d(200,5)
! Compute dv/dy at GYF gridpoints, note remove mean
      do j=2,NY
      do k=0,NZM
      do i=0,NXM
!        S1(i,k,j)=((U2(i,k,j+1)-CR2(0,0,j+1))-(U2(i,k,j)-CR2(0,0,j)))
!     &            /GYF(j)
       S1(i,k,j)=((U2(i,k,j+1)-dble(CR2(0,0,j+1)))
     &        -(U2(i,k,j-1)-dble(CR2(0,0,j-1))))
     &            /(GY(j+1)-GY(j-1))
      end do
      end do
      end do
      do j=2,NY
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+(S1(i,k,j)**2.0)
!        if(k==10) then
        ep2d(i,j)=ep2d(i,j)+(S1(i,k,j)**2.0) 
!        endif
      end do
      end do
      end do
! Compute dw/dy at GYF gridpoints, note remove mean
      do j=2,NY
      do k=0,NZM
      do i=0,NXM
!        S1(i,k,j)=((U3(i,k,j+1)-CR3(0,0,j+1))
!     &      -(U3(i,k,j-1)-CR3(0,0,j-1)))/(GYF(j+1)+GY(j+1))
         S1(i,k,j)=((U3(i,k,j)-dble(CR3(0,0,j)))
     &          -(U3(i,k,j-1)-dble(CR3(0,0,j-1))))
     &             /DY(j)

      end do
      end do
      end do
      do j=2,NY
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+(S1(i,k,j)**2.0)
!        if (k==10) then
         ep2d(i,j)=ep2d(i,j)+(S1(i,k,j)**2.0)
!        endif
      end do
      end do
      end do
      write(66,*) ep2d(200,5)
! Store dv/dz in CF1
      do j=1,NY
      do k=0,TNKZ
      do i=0,NKX
!        CF1(i,k,j)=CIKZ(k)*(CR2(i,k,j)+CR2(i,k,j+1))/2.0
         CF1(i,k,j)=CIKZ(k)*CR2(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xz_to_physical(CF1,F1,0,NY+1)
      do j=2,NY
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+(F1(i,k,j)**2.0)
!        if (k==10) then
         ep2d(i,j)=ep2d(i,j)+(F1(i,k,j)**2.0)
!        endif
      end do
      end do
      end do
      write(66,*) ep2d(200,5)
! Store dw/dz in CS1
      do j=1,NY
      do k=0,TNKZ
      do i=0,NKX
        CS1(i,k,j)=CIKZ(k)*CR3(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=2,NY
      do k=0,NZM
      do i=0,NXM
!        epsilon(j)=epsilon(j)+(S1(i,k,j)**2.0)
         epsilon(j)=epsilon(j)+(DYF(j-1)*S1(i,k,j)**2.d0
     &             +DYF(j)*S1(i,k,j-1)**2.d0)/(2.d0*DY(j))
!         if(k==10) then
         ep2d(i,j)=ep2d(i,j)+(DYF(j-1)*S1(i,k,j)**2.d0
     &             +DYF(j)*S1(i,k,j-1)**2.d0)/(2.d0*DY(j))
!         endif
      end do
      end do
      end do
      write(66,*) ep2d(200,5)
      do j=2,NY
        do i=0,NXM
        ep2d(i,j)=-NU*ep2d(i,j)/float(NZ)
        enddo
      enddo
      do i=1,NX
          ep2d(i,1)=ep2d(i,2)
          ep2d(i,NY+1)=ep2d(i,NY)
      end do
      write(66,*) ep2d(200,1)
      do j=2,NY
        tke_7(j)=-NU*epsilon(j)/float(NX*NZ)
      end do
        tke_7(NY+1)=tke_7(NY)
        tke_7(1)=tke_7(2)
! Get the individual terms of the viscous contribution
! 1,1
      do j=1,NY
        do k=0,TNKZ
          do i=0,NKX
            CS1(i,k,j)=-1.d0*KX2(I)*CR1(i,k,j)
          end do
        end do
      end do
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=S1(i,k,j)*(U1(i,k,j)-dble(CR1(0,0,j)))
          end do
        end do
      end do
      do j=2,NY
        tke_7_1(j)=NU*SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do
      tke_7_1(1)=0.d0
! 1,2
      do j=1,NY+1
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=((U1(i,k,j)-dble(CR1(0,0,j)))
     &               -(U1(i,k,j-1)-dble(CR1(0,0,j-1))))/DY(J)
          end do
        end do
      end do
      do j=1,NY  
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=(U1(i,k,j)-dble(CR1(0,0,j)))
     &               *(S1(i,k,j+1)-S1(i,k,j))/DYF(J)
          end do
        end do
      end do
      do j=2,NY
        tke_7_1(j)=tke_7_1(j)+NU*SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ) 
      end do
! 1,3
      do j=1,NY
        do k=0,TNKZ
          do i=0,NKX
            CS1(i,k,j)=-1.d0*KZ2(K)*CR1(i,k,j)
          end do
        end do
      end do
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=S1(i,k,j)*(U1(i,k,j)-dble(CR1(0,0,j)))
          end do
        end do
      end do
      do j=2,NY
        tke_7_1(j)=tke_7_1(j)+NU*SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do
! 2,1
      do j=2,NY
        do k=0,TNKZ
          do i=0,NKX
            CS1(i,k,j)=-1.d0*KX2(I)*CR2(i,k,j)
          end do
        end do
      end do
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=2,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=S1(i,k,j)*(U2(i,k,j)-dble(CR2(0,0,j)))
          end do
        end do
      end do
      do j=2,NY
        tke_7_2(j)=NU*SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do
! 2,2
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=((U2(i,k,j+1)-dble(CR2(0,0,j+1)))
     &               -(U2(i,k,j)-dble(CR2(0,0,j))))/DYF(J)
          end do
        end do
      end do
      do j=NY,2,-1
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=(U2(i,k,j)-dble(CR2(0,0,j)))
     &               *(S1(i,k,j)-S1(i,k,j-1))/DY(J)
          end do
        end do
      end do
      do j=2,NY
        tke_7_2(j)=tke_7_2(j)+NU*SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do
! 2,3
      do j=2,NY
        do k=0,TNKZ
          do i=0,NKX
            CS1(i,k,j)=-1.d0*KZ2(K)*CR2(i,k,j)
          end do
        end do
      end do
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=2,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=S1(i,k,j)*(U2(i,k,j)-dble(CR2(0,0,j)))
          end do
        end do
      end do
      do j=2,NY
        tke_7_2(j)=tke_7_2(j)+NU*SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do
! 3,1
      do j=1,NY
        do k=0,TNKZ
          do i=0,NKX
            CS1(i,k,j)=-1.d0*KX2(I)*CR3(i,k,j)
          end do
        end do
      end do
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=S1(i,k,j)*(U3(i,k,j)-dble(CR3(0,0,j)))
          end do
        end do
      end do
      do j=2,NY
        tke_7_3(j)=NU*SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do
      tke_7_3(1)=0.d0
! 1,2
      do j=1,NY+1
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=((U3(i,k,j)-dble(CR3(0,0,j)))
     &               -(U3(i,k,j-1)-dble(CR3(0,0,j-1))))/DY(J)
          end do
        end do
      end do
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=(U3(i,k,j)-dble(CR3(0,0,j)))
     &               *(S1(i,k,j+1)-S1(i,k,j))/DYF(J)
          end do
        end do
      end do
      do j=2,NY
        tke_7_3(j)=tke_7_3(j)+NU*SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do
! 1,3
      do j=1,NY
        do k=0,TNKZ
          do i=0,NKX
            CS1(i,k,j)=-1.d0*KZ2(K)*CR3(i,k,j)
          end do
        end do
      end do
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=S1(i,k,j)*(U3(i,k,j)-dble(CR3(0,0,j)))
          end do
        end do
      end do
      do j=2,NY
        tke_7_3(j)=tke_7_3(j)+NU*SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do

! Get the components of the SGS eddy viscosity terms 
!   1,1
      do j=1,NY
        do k=0,TNKZ
          do i=0,NKX
            CS1(i,k,j)=CIKX(I)*CR1(i,k,j)
          end do
        end do
      end do
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=S1(i,k,j)*0.5d0*(NU_T(i,k,j+1)+NU_T(i,k,j))
          end do
        end do
      end do
      call fft_xz_to_fourier(S1,CS1,0,NY+1)
      do j=1,NY
        do k=0,TNKZ
          do i=0,NKX
            CS1(i,k,j)=CIKX(I)*CS1(i,k,j)
          end do
        end do
      end do
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=1,NY
        do k=0,NZM 
          do i=0,NXM
            S1(i,k,j)=S1(i,k,j)*(U1(i,k,j)-dble(CR1(0,0,j)))
          end do 
        end do
      end do
      do j=1,NY
        tke_8_1(j)=SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do 
! 1,2
      do j=1,NY+1
        do k=0,NZM
          do i=0,NXM 
            S1(i,k,j)=NU_T(i,k,j)*((U1(i,k,j)-U1(i,k,j-1))/(DY(J))) 
          end do
        end do
      end do
      do j=1,NY+1
        S1_mean(j)=SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)  
      end do
! Subtract the mean value
      do j=1,NY+1
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=S1(i,k,j)-S1_mean(j)
          end do
        end do
      end do
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=(U1(i,k,j)-dble(CR1(0,0,j)))
     &             *(S1(i,k,j+1)-S1(i,k,j))/DYF(J)
          end do
        end do
      end do
      do j=1,NY
        tke_8_1(j)=tke_8_1(j)+SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do
!   1,3
      do j=1,NY
        do k=0,TNKZ
          do i=0,NKX
            CS1(i,k,j)=CIKZ(K)*CR1(i,k,j)
          end do
        end do
      end do
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=S1(i,k,j)*0.5d0*(NU_T(i,k,j+1)+NU_T(i,k,j))
          end do
        end do
      end do
      call fft_xz_to_fourier(S1,CS1,0,NY+1)
      do j=1,NY
        do k=0,TNKZ
          do i=0,NKX
            CS1(i,k,j)=CIKZ(K)*CS1(i,k,j)
          end do
        end do
      end do
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=S1(i,k,j)*(U1(i,k,j)-dble(CR1(0,0,j)))
          end do
        end do
      end do
      do j=1,NY
        tke_8_1(j)=tke_8_1(j)+SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do
!   2,1
      do j=2,NY
        do k=0,TNKZ
          do i=0,NKX
            CS1(i,k,j)=CIKX(I)*CR2(i,k,j)
          end do
        end do
      end do
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=2,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=S1(i,k,j)*NU_T(i,k,j)
          end do
        end do
      end do
      call fft_xz_to_fourier(S1,CS1,0,NY+1)
      do j=2,NY
        do k=0,TNKZ
          do i=0,NKX
            CS1(i,k,j)=CIKX(I)*CS1(i,k,j)
          end do
        end do
      end do
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=2,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=S1(i,k,j)*(U2(i,k,j)-dble(CR2(0,0,j)))
          end do
        end do
      end do
      do j=2,NY
        tke_8_2(j)=SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do
! 2,2
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=0.5d0*(NU_T(i,k,j+1)+NU_T(i,k,j))
     &              *((U2(i,k,j+1)-U2(i,k,j))/(DYF(J)))
          end do
        end do
      end do
      do j=1,NY
        S1_mean(j)=SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do
! Subtract the mean value
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=S1(i,k,j)-S1_mean(j)
          end do
        end do
      end do
      do j=NY,2,-1
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=(U2(i,k,j)-dble(CR2(0,0,j)))
     &             *(S1(i,k,j)-S1(i,k,j-1))/DY(J)
          end do
        end do
      end do
      do j=2,NY
        tke_8_2(j)=tke_8_2(j)+SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do
!   2,3
      do j=2,NY
        do k=0,TNKZ
          do i=0,NKX
            CS1(i,k,j)=CIKZ(K)*CR2(i,k,j)
          end do
        end do
      end do
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=2,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=S1(i,k,j)*NU_T(i,k,j)
          end do
        end do
      end do
      call fft_xz_to_fourier(S1,CS1,0,NY+1)
      do j=2,NY
        do k=0,TNKZ
          do i=0,NKX
            CS1(i,k,j)=CIKZ(K)*CS1(i,k,j)
          end do
        end do
      end do
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=2,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=S1(i,k,j)*(U2(i,k,j)-dble(CR2(0,0,j)))
          end do
        end do
      end do
      do j=2,NY
        tke_8_2(j)=tke_8_2(j)+SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do
!   3,1
      do j=1,NY
        do k=0,TNKZ
          do i=0,NKX
            CS1(i,k,j)=CIKX(I)*CR3(i,k,j)
          end do
        end do
      end do
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=S1(i,k,j)*0.5d0*(NU_T(i,k,j+1)+NU_T(i,k,j))
          end do
        end do
      end do
      call fft_xz_to_fourier(S1,CS1,0,NY+1)
      do j=1,NY
        do k=0,TNKZ
          do i=0,NKX
            CS1(i,k,j)=CIKX(I)*CS1(i,k,j)
          end do
        end do
      end do
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=S1(i,k,j)*(U3(i,k,j)-dble(CR3(0,0,j)))
          end do
        end do
      end do
      do j=1,NY
        tke_8_3(j)=SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do
! 3,2
      do j=1,NY+1
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=NU_T(i,k,j)*((U3(i,k,j)-U3(i,k,j-1))/(DY(J)))
          end do
        end do
      end do
      do j=1,NY+1
        S1_mean(j)=SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do
! Subtract the mean value
      do j=1,NY+1
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=S1(i,k,j)-S1_mean(j)
          end do
        end do
      end do
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=(U3(i,k,j)-dble(CR3(0,0,j)))
     &             *(S1(i,k,j+1)-S1(i,k,j))/DYF(J)
          end do
        end do
      end do
      do j=1,NY
        tke_8_3(j)=tke_8_3(j)+SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do
!   3,3
      do j=1,NY
        do k=0,TNKZ
          do i=0,NKX
            CS1(i,k,j)=CIKZ(K)*CR3(i,k,j)
          end do
        end do
      end do
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=S1(i,k,j)*0.5d0*(NU_T(i,k,j+1)+NU_T(i,k,j))
          end do
        end do
      end do
      call fft_xz_to_fourier(S1,CS1,0,NY+1)
      do j=1,NY
        do k=0,TNKZ
          do i=0,NKX
            CS1(i,k,j)=CIKZ(K)*CS1(i,k,j)
          end do
        end do
      end do
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            S1(i,k,j)=S1(i,k,j)*(U3(i,k,j)-dble(CR3(0,0,j)))
          end do
        end do
      end do
      do j=1,NY
        tke_8_3(j)=tke_8_3(j)+SUM(S1(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do



! Write out the mean statistics at each time
c      open(45,file='tke.txt',form='formatted',status='unknown')

      open(45,file=FNAME,form='formatted',status='unknown',
     + position='append')


      write(45,*) TIME_STEP,TIME,DELTA_T
      do j=1,NY
        write(45,401) j,GY(J),tke_1(j),tke_2(j),tke_3(j),tke_4(j)
     &          ,(tke_5(j,n),n=1,N_TH),tke_6_1(j)
     &       ,tke_6_1_1(j),tke_6_1_2(j),tke_6_1_3(j)
     &          ,tke_6_2(j)
     &       ,tke_6_2_1(j),tke_6_2_2(j),tke_6_2_3(j)
     &          ,tke_7(j)
     &          ,tke_3_1(j),tke_3_2(j),tke_3_3(j)
     &          ,tke_7_1(j),tke_7_2(j),tke_7_3(j)
     &          ,tke_8_1(j),tke_8_2(j),tke_8_3(j)
      end do
      close(45)
    
401   format(I3,' ',24(F20.9,' '))

      return 
      end


c----*|--------------------------------------------
      SUBROUTINE SAVE_PLANE_VEL
c----*|---------------------------------------------
!
!     subroutine to save planes
!      
      INCLUDE 'header'
      INCLUDE 'mpif.h'
      INCLUDE 'header_mpi'

      INTEGER i,j,k
      CHARACTER*38 filn_IW_1,filn_IW_2,
     &             filn_IW_xz_1,filn_IW_xz_2,
     &             filn_IW_xz_3

      CHARACTER*38 filn_vel_xy,filn_vel_yz 
      CHARACTER*38 filn_vel_xz_1,filn_vel_xz_2,
     &             filn_vel_xz_3

      i = rank + 1

      j = TIME_STEP

      filn_IW_1 = 'plane_data/data_IW_xy_1'
     &         //CHAR(MOD(I,100)/10+48)
     &        //CHAR(MOD(I,10)+48) // '_'
     &        //CHAR(MOD(j,10000000)/1000000+48)
     &        //CHAR(MOD(j,1000000)/100000+48)
     &        //CHAR(MOD(j,100000)/10000+48)
     &        //CHAR(MOD(j,10000)/1000+48)
     &        //CHAR(MOD(j,1000)/100+48)
     &        //CHAR(MOD(j,100)/10+48)
     &        //CHAR(MOD(j,10)+48) //
     &        '.pln'

      filn_IW_2 = 'plane_data/data_IW_xy_2'
     &         //CHAR(MOD(I,100)/10+48)
     &        //CHAR(MOD(I,10)+48) // '_'
     &        //CHAR(MOD(j,10000000)/1000000+48)
     &        //CHAR(MOD(j,1000000)/100000+48)
     &        //CHAR(MOD(j,100000)/10000+48)
     &        //CHAR(MOD(j,10000)/1000+48)
     &        //CHAR(MOD(j,1000)/100+48)
     &        //CHAR(MOD(j,100)/10+48)
     &        //CHAR(MOD(j,10)+48) //
     &        '.pln'



      filn_vel_xy = 'plane_data/data_vel_xy'
     &         //CHAR(MOD(I,100)/10+48)
     &        //CHAR(MOD(I,10)+48) // '_'
     &        //CHAR(MOD(j,10000000)/1000000+48)
     &        //CHAR(MOD(j,1000000)/100000+48)
     &        //CHAR(MOD(j,100000)/10000+48)
     &        //CHAR(MOD(j,10000)/1000+48)
     &        //CHAR(MOD(j,1000)/100+48)
     &        //CHAR(MOD(j,100)/10+48)
     &        //CHAR(MOD(j,10)+48) //
     &        '.pln'

       filn_vel_yz = 'plane_data/data_vel_yz'
     &         //CHAR(MOD(I,100)/10+48)
     &        //CHAR(MOD(I,10)+48) // '_'
     &        //CHAR(MOD(j,10000000)/1000000+48)
     &        //CHAR(MOD(j,1000000)/100000+48)
     &        //CHAR(MOD(j,100000)/10000+48)
     &        //CHAR(MOD(j,10000)/1000+48)
     &        //CHAR(MOD(j,1000)/100+48)
     &        //CHAR(MOD(j,100)/10+48)
     &        //CHAR(MOD(j,10)+48) //
     &        '.pln'


      open(443,file=filn_IW_1,status='unknown',form='unformatted')
      open(444,file=filn_IW_2,status='unknown',form='unformatted')

      open(445,file=filn_vel_xy,status='unknown',form='unformatted')
      open(449,file=filn_vel_yz,status='unknown',form='unformatted')

       

           write(443)TIME
           write(443) (dble(CR1(0,0,j)),j=1,NY), (((U2(i,0,j)-
     &      dble(CR2(0,0,j))),I=0,NXM),J=1,NY),(((((U2(i,0,j+1)
     &     -dble(CR2(0,0,j+1)))
     &    -(U2(i,0,j)-dble(CR2(0,0,j))))/(DYF(j))),I=0,NXM),J=1,NY)

           close(443)

           write(444)TIME
           write(444) (dble(CR1(0,0,j)),j=1,NY), (((U2(i,NZM,j)-
     &      dble(CR2(0,0,j))),I=0,NXM),J=1,NY),(((((U2(i,NZM,j+1)
     &     -dble(CR2(0,0,j+1)))
     &    -(U2(i,NZM,j)-dble(CR2(0,0,j))))/(DYF(j))),I=0,NXM),J=1,NY)

           close(444)

           write(445)TIME
           write(445)((P(i,0,j),i=0,NXM),j=1,NY),
     &     ((U1(i,0,j),i=0,NXM),j=1,NY),
     &     ((U3(i,0,j),i=0,NXM),j=1,NY),
     &     ((S1(i,0,j),i=0,NXM),j=1,NY),
     &     ((P(i,NZM,j),i=0,NXM),j=1,NY),
     &     ((U1(i,NZM,j),i=0,NXM),j=1,NY),
     &     ((U3(i,NZM,j),i=0,NXM),j=1,NY),
     &     ((S1(i,NZM,j),i=0,NXM),j=1,NY) 
  
           close(445)
         

           write(449)TIME
           write(449)((P(0,k,j),k=0,NZM),j=1,NY),
     &     ((U1(0,k,j),k=0,NZM),j=1,NY),
     &     ((U2(0,k,j),k=0,NZM),j=1,NY),
     &     ((U3(0,k,j),k=0,NZM),j=1,NY),
     &     ((P(NXM,k,j),k=0,NZM),j=1,NY),
     &     ((U1(NXM,k,j),k=0,NZM),j=1,NY),
     &     ((U1(NXM,k,j),k=0,NZM),j=1,NY),
     &     ((U3(NXM,k,j),k=0,NZM),j=1,NY)

           close(449)

           IF    (RANK .eq. 0) THEN   

          i = RANK + 1
          j = TIME_STEP      

           filn_vel_xz_1 = 'plane_data/data_vel_xz_1'
     &         //CHAR(MOD(I,100)/10+48)
     &        //CHAR(MOD(I,10)+48) // '_'
     &        //CHAR(MOD(j,10000000)/1000000+48)
     &        //CHAR(MOD(j,1000000)/100000+48)
     &        //CHAR(MOD(j,100000)/10000+48)
     &        //CHAR(MOD(j,10000)/1000+48)
     &        //CHAR(MOD(j,1000)/100+48)
     &        //CHAR(MOD(j,100)/10+48)
     &        //CHAR(MOD(j,10)+48) //
     &        '.pln'

           filn_IW_xz_1 = 'plane_data/data_IW_xz_1'
     &         //CHAR(MOD(I,100)/10+48)
     &        //CHAR(MOD(I,10)+48) // '_'
     &        //CHAR(MOD(j,10000000)/1000000+48)
     &        //CHAR(MOD(j,1000000)/100000+48)
     &        //CHAR(MOD(j,100000)/10000+48)
     &        //CHAR(MOD(j,10000)/1000+48)
     &        //CHAR(MOD(j,1000)/100+48)
     &        //CHAR(MOD(j,100)/10+48)
     &        //CHAR(MOD(j,10)+48) //
     &        '.pln'

         
        j =  80
        open(446,file=filn_vel_xz_1,status='unknown',form='unformatted')

        open(646,file=filn_IW_xz_1,status='unknown',form='unformatted')

         write(446)TIME
         write(446)((P(i,k,j),i=0,NXM),k=0,NZM),
     &    ((U1(i,k,j),i=0,NXM),k=0,NZM),
     &    ((U2(i,k,j),i=0,NXM),k=0,NZM),
     &    ((U3(i,k,j),i=0,NXM),k=0,NZM)
    
         close(446)

         write(646)TIME
         write(646) dble(CR1(0,0,j)), (((U2(i,k,j)-
     &      dble(CR2(0,0,j))),I=0,NXM),k=0,NZM),(((((U2(i,k,j+1)
     &     -dble(CR2(0,0,j+1)))
     &    -(U2(i,k,j)-dble(CR2(0,0,j))))/(DYF(j))),I=0,NXM),k=0,NZM)

          close(646)  

           ELSEIF (RANK .eq. 2 )THEN

         i = RANK + 1
          j = TIME_STEP

           filn_vel_xz_2 = 'plane_data/data_vel_xz_2'
     &         //CHAR(MOD(I,100)/10+48)
     &        //CHAR(MOD(I,10)+48) // '_'
     &        //CHAR(MOD(j,10000000)/1000000+48)
     &        //CHAR(MOD(j,1000000)/100000+48)
     &        //CHAR(MOD(j,100000)/10000+48)
     &        //CHAR(MOD(j,10000)/1000+48)
     &        //CHAR(MOD(j,1000)/100+48)
     &        //CHAR(MOD(j,100)/10+48)
     &        //CHAR(MOD(j,10)+48) //
     &        '.pln'

           filn_IW_xz_2 = 'plane_data/data_IW_xz_2'
     &         //CHAR(MOD(I,100)/10+48)
     &        //CHAR(MOD(I,10)+48) // '_'
     &        //CHAR(MOD(j,10000000)/1000000+48)
     &        //CHAR(MOD(j,1000000)/100000+48)
     &        //CHAR(MOD(j,100000)/10000+48)
     &        //CHAR(MOD(j,10000)/1000+48)
     &        //CHAR(MOD(j,1000)/100+48)
     &        //CHAR(MOD(j,100)/10+48)
     &        //CHAR(MOD(j,10)+48) //
     &        '.pln' 
           

        j = 17
        open(447,file=filn_vel_xz_2,status='unknown',form='unformatted')

        open(647,file=filn_IW_xz_2,status='unknown',form='unformatted')

         write(447)TIME
         write(447)((P(i,k,j),i=0,NXM),k=0,NZM),
     &    ((U1(i,k,j),i=0,NXM),k=0,NZM),
     &    ((U2(i,k,j),i=0,NXM),k=0,NZM),
     &    ((U3(i,k,j),i=0,NXM),k=0,NZM)
   
         close(447)

         write(647)TIME
         write(647) dble(CR1(0,0,j)), (((U2(i,k,j)-
     &      dble(CR2(0,0,j))),I=0,NXM),k=0,NZM),(((((U2(i,k,j+1)
     &     -dble(CR2(0,0,j+1)))
     &    -(U2(i,k,j)-dble(CR2(0,0,j))))/(DYF(j))),I=0,NXM),k=0,NZM)

          close(647)

           ELSEIF (RANK .eq. 3 )THEN
        
         i = RANK + 1
         j = TIME_STEP

           filn_vel_xz_3 = 'plane_data/data_vel_xz_3'
     &         //CHAR(MOD(I,100)/10+48)
     &        //CHAR(MOD(I,10)+48) // '_'
     &        //CHAR(MOD(j,10000000)/1000000+48)
     &        //CHAR(MOD(j,1000000)/100000+48)
     &        //CHAR(MOD(j,100000)/10000+48)
     &        //CHAR(MOD(j,10000)/1000+48)
     &        //CHAR(MOD(j,1000)/100+48)
     &        //CHAR(MOD(j,100)/10+48)
     &        //CHAR(MOD(j,10)+48) //
     &        '.pln'

            filn_IW_xz_3 = 'plane_data/data_IW_xz_3'
     &         //CHAR(MOD(I,100)/10+48)
     &        //CHAR(MOD(I,10)+48) // '_'
     &        //CHAR(MOD(j,10000000)/1000000+48)
     &        //CHAR(MOD(j,1000000)/100000+48)
     &        //CHAR(MOD(j,100000)/10000+48)
     &        //CHAR(MOD(j,10000)/1000+48)
     &        //CHAR(MOD(j,1000)/100+48)
     &        //CHAR(MOD(j,100)/10+48)
     &        //CHAR(MOD(j,10)+48) //
     &        '.pln'

        j = 27

        open(448,file=filn_vel_xz_3,status='unknown',form='unformatted')

        open(648,file=filn_IW_xz_3,status='unknown',form='unformatted')

         write(448)TIME
         write(448)((P(i,k,j),i=0,NXM),k=0,NZM),
     &    ((U1(i,k,j),i=0,NXM),k=0,NZM),
     &    ((U2(i,k,j),i=0,NXM),k=0,NZM),
     &    ((U3(i,k,j),i=0,NXM),k=0,NZM)

         close(448)


         write(648)TIME
         write(648) dble(CR1(0,0,j)), (((U2(i,k,j)-
     &      dble(CR2(0,0,j))),I=0,NXM),k=0,NZM),(((((U2(i,k,j+1)
     &     -dble(CR2(0,0,j+1)))
     &    -(U2(i,k,j)-dble(CR2(0,0,j))))/(DYF(j))),I=0,NXM),k=0,NZM)

          close(648)

 
           ENDIF  

          return
          end
 




      SUBROUTINE SAVE_PLANE
!      
      INCLUDE 'header'
      INCLUDE 'mpif.h'
      INCLUDE 'header_mpi'

      INTEGER i, k

!!!!!!!  Initialization !!!!!!!!!!!!!!!!
          U1_t(:,:) = 0.0d0 ;
          U2_t(:,:) = 0.0d0 ;
          U3_t(:,:) = 0.0d0 ;
          TH_t(:,:) = 0.0d0 ;
          S_t(:,:)  = 0.0d0 ;
         tket(:,:) = 0.0d0 ;
          disst(:,:)=0.0d0;
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          U1_ty(:,:) = 0.0d0 ;
          U2_ty(:,:) = 0.0d0 ;
          U3_ty(:,:) = 0.0d0 ;
          TH_ty(:,:) = 0.0d0 ;
          S_ty(:,:)  = 0.0d0 ;


       do i=1,NX
         CALL MPI_INST_PROF(U1(i,0,1:NY),U1_t(1:NY_TOT,i))
         CALL MPI_INST_PROF(U2(i,0,1:NY),U2_t(1:NY_TOT,i))
         CALL MPI_INST_PROF(U3(i,0,1:NY),U3_t(1:NY_TOT,i))
         CALL MPI_INST_PROF(TH(i,0,1:NY,1),TH_t(1:NY_TOT,i))
         CALL MPI_INST_PROF(TH(i,0,1:NY,2),S_t(1:NY_TOT,i))
         CALL MPI_INST_PROF(tke2d(i,1:NY),tket(1:NY_TOT,i))
         CALL MPI_INST_PROF(ep2d(i,1:NY),disst(1:NY_TOT,i))
         CALL MPI_INST_PROF(bp2d(i,1:NY),buot(1:NY_TOT,i))
       enddo

       do k=1,NZ
         CALL MPI_INST_PROF(U1(0,k,1:NY),U1_ty(1:NY_TOT,k))
         CALL MPI_INST_PROF(U2(0,k,1:NY),U2_ty(1:NY_TOT,k))
         CALL MPI_INST_PROF(U3(0,k,1:NY),U3_ty(1:NY_TOT,k))
         CALL MPI_INST_PROF(TH(0,k,1:NY,1),TH_ty(1:NY_TOT,k))
         CALL MPI_INST_PROF(TH(0,k,1:NY,2),S_ty(1:NY_TOT,k))
!        CALL MPI_INST_PROF(tke2d(k,1:NY),tke_ty(1:NY_TOT,k))
       enddo
      
       write(6,*) 'dissipation is', disst(400,5),ep2d(20,5)
     

       do i=1,NX
         do k=1,NZ
          U2m(k,i)=((TH(i,k,2,1)-TH(i,k,1,1))/(GYF(2)-GYF(1)))*
     &   (C_sp_heat/L_heat)*(NU/Pr(1))*(1.0d0/0.92)
          TH_f(k,i)=(TH(i,k,2,1)-TH(i,k,1,1))
          S_f(k,i) =(TH(i,k,2,2)-TH(i,k,1,2))
         enddo
       enddo

       IF (RANK == 0) THEN

        CALL plane_parav_vel_wr_xz(TIME_STEP)
        CALL plane_parav_vel_wr_yx(TIME_STEP)
        CALL plane_parav_vel_wr_yz(TIME_STEP)
!        CALL plane_parav_tke_xz(TIME_STEP)
       ENDIF

       return
       end
      SUBROUTINE SAVE_PLANE_3D
c----*|---------------------------------------------
!
!     subroutine to save planes
!      
      INCLUDE 'header_les'
      INCLUDE 'mpif.h'
      INCLUDE 'header_mpi'
      

      INTEGER i,j,k
      CHARACTER*38 filn_3D

      i = rank + 1

      j = TIME_STEP

      filn_3D = 'plane_3D/data_3D_'
     &         //CHAR(MOD(I,100)/10+48)
     &        //CHAR(MOD(I,10)+48) // '_'
     &        //CHAR(MOD(j,10000000)/1000000+48)
     &        //CHAR(MOD(j,1000000)/100000+48)
     &        //CHAR(MOD(j,100000)/10000+48)
     &        //CHAR(MOD(j,10000)/1000+48)
     &        //CHAR(MOD(j,1000)/100+48)
     &        //CHAR(MOD(j,100)/10+48)
     &        //CHAR(MOD(j,10)+48) //
     &        '.pln'


      open(543,file=filn_3D,status='unknown',form='unformatted')

 

           write(543)TIME
           write(543)(C_DYN(j),j=1,NY),
     &      (((U1(i,k,j),i=0,NXM),k=0,NZM),j=1,NY), 
     &      (((U2(i,k,j),i=0,NXM),k=0,NZM),j=1,NY),
     &      (((U3(i,k,j),i=0,NXM),k=0,NZM),j=1,NY), 
     &      (((TH(i,k,j,1),i=0,NXM),k=0,NZM),j=1,NY)
             
           close(543)
         
          return
          end 
!modification by mainak for writting into netcdf format
!      SUBROUTINE nc_write


!      INCLUDE 'header'
!      INCLUDE "mpif.h"
!      INCLUDE 'header_mpi'
!      include 'netcdf.inc'


!          INTEGER :: dimid,stat
!          INTEGER :: x_varid, y_varid, varid(11),disid
!          stat=nf_open("profile.nc",NF_WRITE,ncid)
!          if(stat.eq.nf_noerr)then
!            goto 114
!          else
!          endif
!          CALL check(nf_create("profile.nc", NF_CLOBBER,ncid))
!          write(6,*) 'netcdf writing test'
!          CALL check(nf_def_dim(ncid,"Wall_normal",NY*NPROCS,dimid))
!         CALL check(nf_def_var(ncid,"distance",NF_REAL,dimid,disid))
!         CALL check(nf_def_var(ncid,"bf_flo",NF_REAL,dimid,varid(1)))
!          CALL check(nf_def_var(ncid,"temp",NF_REAL,dimid,varid(2)))
!          CALL check(nf_def_var(ncid,"sal",NF_REAL,dimid,varid(3)))
!          CALL check(nf_enddef(ncid))
!114       CALL check(nf_put_var(ncid,disid,GYF_TOT))
!          CALL check(nf_put_var(ncid,varid(1),UME_TOT))
!          CALL check(nf_put_var(ncid,varid(2),TH_TOT))
!          CALL check(nf_put_var(ncid,varid(3),SAL_TOT))
!          CALL check(nf_close(ncid))
!      END SUBROUTINE nc_write
          
!       subroutine check(status)
!       include 'netcdf.inc'
!       integer, intent ( in) :: status

!        if(status /= nf_noerr) then
!          write(6,*) 'somethings wrong'
!          print *, trim(nf_strerror(status))
!          stop "Stopped"
!        end if
       ! end subroutine check

      subroutine plane_parav_vel_wr_xz(jj)

      
       INCLUDE 'header'

       integer i,j,k,jj
       real*8     zcor(NX)
       real*4     var_1(1:3*NX*NY_TOT)
       real*8     salt(1:NY_TOT,1:NX)
       !FILENAMES
       character(len=25) :: ss,basename

       integer :: s11,np2,np1
       CHARACTER*38 fil_tim
       np1=NX
       np2= NY_TOT
       do i= 1,np1
         zcor(i)=LX/dble(NX-1)*dble(i-1)
       enddo


        fil_tim   = 'plane_data/data_xz_'
     &        //CHAR(MOD(jj,10000000)/1000000+48)
     &        //CHAR(MOD(jj,1000000)/100000+48)
     &        //CHAR(MOD(jj,100000)/10000+48)
     &        //CHAR(MOD(jj,10000)/1000+48)
     &        //CHAR(MOD(jj,1000)/100+48)
     &        //CHAR(MOD(jj,100)/10+48)
     &        //CHAR(MOD(jj,10)+48) //
     &        '.vtk'

       var_1(:) = 0.0d0 ;

         k = 1
         do i = 1,np1
         do j = 1,np2
         var_1(k) =  real(U2_t(j,i))
         k = k+1
         var_1(k) =  real(U3_t(j,i))
         k = k+1
         var_1(k) =  real(U1_t(j,i))
         k = k+1
        enddo
        enddo
        
         do i=1,np1
          do j=1,np2
            salt(j,i)=Sal_0-Grad_rho(2)*(Cos(ANG_BETA)*(GYF_TOT(np2)
     & -GYF_TOT(j))+Sin(ANG_BETA)*(zcor(i)-zcor(1)))+S_T(j,i)
          enddo
         enddo

      open(unit=13,file=fil_tim,access='stream',  form='unformatted'
     & ,status='unknown', convert='big_endian',iostat=s11)



!HEADER: note termination with char(10)
        
        write(13) "# vtk DataFile Version 3.0"//char(10)
        write(13) trim(BaseName)//char(10)
        write(13) "BINARY"//char(10)
        write(13) "DATASET RECTILINEAR_GRID"//char(10)

        write(ss,fmt='(A10,3I5)') "DIMENSIONS",np2,1,np1
        write(13) ss//char(10)
!        write(ss,fmt='(A4,I6,A6)') "TIME",1," float"
!        write(13) char(10)//ss//char(10)
!        write(13) real(TIME)
        
        write(ss,fmt='(A13,I6,A6)') "X_COORDINATES",np2," float"
        write(13) char(10)//ss//char(10)
        do i = 1,np2
        write(13)real(GYF_TOT(i) )
        enddo
       
        write(ss,fmt='(A13,I6,A6)') "Y_COORDINATES",1," float"
        write(13) char(10)//ss//char(10)
        write(13) 0.01
        
        write(ss,fmt='(A13,I6,A6)') "Z_COORDINATES",np1," float"
        write(13) char(10)//ss//char(10)
        do j = 1,np1
        write(13) real(zcor(j))
        enddo
        
        write(ss,fmt='(A10,I15)') "POINT_DATA",np1*np2
        write(13) char(10)//ss//char(10)
        write(13) "VECTORS velocity_vectors float 1"//char(10)
        write(13) var_1
        write(13) "SCALARS temperature float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(TH_t(1:np2,1:np1))
        write(13) "SCALARS salinity_anom float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(S_t(1:np2,1:np1))
         write(13) "SCALARS salinity float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(salt(1:np2,1:np1))
        close(13)
        return
        end
       
       subroutine plane_parav_vel_wr_yx(jj)

       INCLUDE 'header'

       integer i,j,k,jj
       real*8     ycor(NZ)
       !FILENAMES
       character(len=25) :: ss,basename
       real*4     var_2(1:3*NZ*NY_TOT)
       integer ::  s11,np2,np3
       CHARACTER*38 fil_time
       np3=NZ
       np2= NY_TOT
       do k= 1,np3
         ycor(k)=LZ/dble(NZ-1)*dble(k-1)
       enddo
       var_2(:)=0.0d0 ;
       k = 1
         do i = 1,np3
         do j = 1,np2
         var_2(k) =  real(U2_ty(j,i))
         k = k+1
         var_2(k) =  real(U3_ty(j,i))
         k = k+1
         var_2(k) =  real(U1_ty(j,i))
         k = k+1
        enddo
        enddo
       
       fil_time  = 'plane_data/data_yx'
     &        //CHAR(MOD(jj,10000000)/1000000+48)
     &        //CHAR(MOD(jj,1000000)/100000+48)
     &        //CHAR(MOD(jj,100000)/10000+48)
     &        //CHAR(MOD(jj,10000)/1000+48)
     &        //CHAR(MOD(jj,1000)/100+48)
     &        //CHAR(MOD(jj,100)/10+48)
     &        //CHAR(MOD(jj,10)+48) //
     &        '.vtk'
     
        open(unit=13,file=fil_time,access='stream',  form='unformatted'
     &  ,status='unknown', convert='big_endian',iostat=s11)
    
        write(13) "# vtk DataFile Version 3.0"//char(10)
        write(13) trim(BaseName)//char(10)
        write(13) "BINARY"//char(10)
        write(13) "DATASET RECTILINEAR_GRID"//char(10)

        write(ss,fmt='(A10,3I5)') "DIMENSIONS",np2,np3,1
        write(13) ss//char(10)

        write(ss,fmt='(A13,I6,A6)') "X_COORDINATES",np2," float"
        write(13) char(10)//ss//char(10)
        do i = 1,np2
        write(13)real(GYF_TOT(i) )
        enddo
        
        write(ss,fmt='(A13,I6,A6)') "Y_COORDINATES",np3," float"
        write(13) char(10)//ss//char(10)
        do k = 1,np3
          write(13)real(ycor(k))
        enddo
        
        write(ss,fmt='(A13,I6,A6)') "Z_COORDINATES",1," float"
        write(13) char(10)//ss//char(10)
        write(13) 0.0

        write(ss,fmt='(A10,I15)') "POINT_DATA",np2*np3
        write(13) char(10)//ss//char(10) 
        write(13) "VECTORS velocity_vectors float 1"//char(10)
        write(13) var_2
        write(13) "SCALARS Spanwise_vel float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(U3_ty(1:np2,1:np3))
        write(13) "SCALARS temperature float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(TH_ty(1:np2,1:np3))
        write(13) "SCALARS salinity float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(S_ty(1:np2,1:np3))
        close(13)
        return
        end
 
       subroutine plane_parav_vel_wr_yz(jj)

       INCLUDE 'header'

       integer i,j,k,jj
       real*8     zcor(NX),ycor(NZ)
       integer :: s11,np1,np3
       CHARACTER*38 fil_tme
       character(len=25) :: ss,basename
       np3=NZ
       np1=NX
       do k= 1,np3
         ycor(k)=LZ/dble(NZ-1)*dble(k-1)
       enddo
       do i= 1,np1
         zcor(i)=LX/dble(NX-1)*dble(i-1)
       enddo
        fil_tme  = 'plane_data/data_yz'
     &        //CHAR(MOD(jj,10000000)/1000000+48)
     &        //CHAR(MOD(jj,1000000)/100000+48)
     &        //CHAR(MOD(jj,100000)/10000+48)
     &        //CHAR(MOD(jj,10000)/1000+48)
     &        //CHAR(MOD(jj,1000)/100+48)
     &        //CHAR(MOD(jj,100)/10+48)
     &        //CHAR(MOD(jj,10)+48) //
     &        '.vtk'
     
        open(unit=13,file=fil_tme,access='stream',  form='unformatted'
     &  ,status='unknown', convert='big_endian',iostat=s11)

        write(13) "# vtk DataFile Version 3.0"//char(10)
        write(13) trim(BaseName)//char(10)
        write(13) "BINARY"//char(10)
        write(13) "DATASET RECTILINEAR_GRID"//char(10)

        write(ss,fmt='(A10,3I5)') "DIMENSIONS",1,np3,np1
        write(13) ss//char(10)
        write(ss,fmt='(A13,I6,A6)') "X_COORDINATES",1," float"
        write(13) char(10)//ss//char(10)
        write(13) 0.0108
        write(ss,fmt='(A13,I6,A6)') "Y_COORDINATES",np3," float"
        write(13) char(10)//ss//char(10)
        do k = 1,np3
          write(13)real(ycor(k))
        enddo
        write(ss,fmt='(A13,I6,A6)') "Z_COORDINATES",np1," float"
        write(13) char(10)//ss//char(10)
        do i= 1,np1
        write(13) real(zcor(i))
        enddo
        write(ss,fmt='(A10,I15)') "POINT_DATA",np1*np3
        write(13) char(10)//ss//char(10)
        write(13) "SCALARS Meltrate float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(U2m(1:np3,1:np1))
        write(13) "SCALARS Temp_flux float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(TH_f(1:np3,1:np1))
        write(13) "SCALARS Sal_flux float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(S_f(1:np3,1:np1))
        close(13)
        return
        end
         
          subroutine plane_parav_tke_xz(jj)


        INCLUDE 'header'

        integer i,j,k,jj
        real*8     zcor(NX)
        real*4     var_1(1:3*NX*NY_TOT)
        !FILENAMES
        character(len=25) :: ss,basename

        integer :: s11,np2,np1
        CHARACTER*38 fil_tim
        np1=NX
        np2= NY_TOT
        do i= 1,np1
          zcor(i)=LX/dble(NX-1)*dble(i-1)
        enddo


         fil_tim   = 'plane_data/data_tke_xz_'
     &        //CHAR(MOD(jj,10000000)/1000000+48)
     &        //CHAR(MOD(jj,1000000)/100000+48)
     &        //CHAR(MOD(jj,100000)/10000+48)
     &        //CHAR(MOD(jj,10000)/1000+48)
     &        //CHAR(MOD(jj,1000)/100+48)
     &        //CHAR(MOD(jj,100)/10+48)
     &        //CHAR(MOD(jj,10)+48) //
     &        '.vtk'
     
        open(unit=13,file=fil_tim,access='stream',  form='unformatted'
     & ,status='unknown', convert='big_endian',iostat=s11)



!HEADER: note termination with char(10)

        write(13) "# vtk DataFile Version 3.0"//char(10)
        write(13) trim(BaseName)//char(10)
        write(13) "BINARY"//char(10)
        write(13) "DATASET RECTILINEAR_GRID"//char(10)

        write(ss,fmt='(A10,3I5)') "DIMENSIONS",np2,1,np1
        write(13) ss//char(10)

        write(ss,fmt='(A13,I6,A6)') "X_COORDINATES",np2," float"
        write(13) char(10)//ss//char(10)
        do i = 1,np2
        write(13)real(GYF_TOT(i) )
        enddo

        write(ss,fmt='(A13,I6,A6)') "Y_COORDINATES",1," float"
        write(13) char(10)//ss//char(10)
        write(13) 0.01

        write(ss,fmt='(A13,I6,A6)') "Z_COORDINATES",np1," float"
        write(13) char(10)//ss//char(10)
        do j = 1,np1
        write(13) real(zcor(j))
        enddo

        write(ss,fmt='(A10,I15)') "POINT_DATA",np1*np2
        write(13) char(10)//ss//char(10) 
        write(13) "SCALARS tke float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(tket(1:np2,1:np1))
        write(13) "SCALARS dissip float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(disst(1:np2,1:np1))
        write(13) "SCALARS bu_prod float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(buot(1:np2,1:np1))

        close(13)
        return
        end
