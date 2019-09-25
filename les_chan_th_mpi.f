      subroutine les_chan_th(n)
C This subroutine models the subgridscale terms 
c in the scalar advection equation for scalar number n
C if the computation is to be treated as an LES not a DNS
C This subroutine should be called when the velocity is in fourier space 
C   in the periodic directions
C S1 should contain |S| which was calculated in les_chan

      include 'header_les'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'

      CHARACTER*35 FNAME_TH
      integer i,j,k,l,m,ij,n

      real*8 S1_mean(1:NY)
 

! Variables for Dynamic Smagorinsky model:
      real*8 C_SMAG
      parameter (C_SMAG=0.13d0)
      real*8 DELTA_Y(0:NY+1)
      real*8 alpha,beta 
      real*8 denominator_sum(0:NY+1), numerator_sum(0:NY+1)

! Here, alpha is the test/LES filter width ratio
      parameter (alpha=2.44950)
! beta is the LES/grid filter width ratio
      parameter (beta=1.d0)


      I = RANK+1

           FNAME_TH ='results/mean_les_th_'
     &        //CHAR(MOD(I,100)/10+48)
     &        //CHAR(MOD(I,10)+48) // '.txt'
! First, for all models, apply boundary conditions to the velocity field
! (fill ghost cells) to ensure accurate calculation of gradients
      IF (USE_MPI) THEN
        J1=JSTART
        J2=JEND 
        CALL APPLY_BC_VEL_MPI
      ELSE
        CALL APPLY_BC_VEL_LOWER
        CALL APPLY_BC_VEL_UPPER
      END IF

      if (LES_MODEL_TYPE.EQ.1) then
!     Constant Smagorinsky model

! First, compute the rate of strain tensor S_ij

      call compute_scalar_grad(n)

! Now, compute |S|*dTH/dx_i, storing in Sij
! First compute at GYF points 
      DO J=0,NY+1
        DO K=0,NZM
          DO I=0,NXM
            Sij(I,K,J,1)=S1(I,K,J)*Sij(I,K,J,1)
            Sij(I,K,J,3)=S1(I,K,J)*Sij(I,K,J,3)
          END DO
        END DO
      END DO
! Convert |S|*S_ij to Fourier space
      ij=1 
      CALL FFT_XZ_TO_FOURIER(Sij(0,0,0,ij),CSij(0,0,0,ij),0,NY+1)
      ij=3
      CALL FFT_XZ_TO_FOURIER(Sij(0,0,0,ij),CSij(0,0,0,ij),0,NY+1)

! Sij(:,:,:,2) is added through an implicit eddy viscosity
      DO J=1,NY+1
        DO K=0,TNKZ
          DO I=0,NKX
             CSij(I,K,J,2)=0.d0
          END DO
        END DO
      END DO

! We now have |S|*dTH/dx_i stored in Sij(:,:,:,1..3) in Physical space

! Compute the filter lengthscale
! Absorb -2.d0*C_SMAG**2.d0 here for effienciency
      DO J=1,NY
! At GYF points:
! Constant Smagorinsky
        DELTA_YF(J)=-C_SMAG**2.d0
     &     *(DX(1)*beta*DYF(J)*2.d0*DZ(1)*beta)**(2.d0/3.d0)
! Wall Damping
!        DELTA_YF(J)=
!     &    -1.d0*(0.1d0*(1.d0-exp((-GYF(J)/(NU*25.d0))**3.d0)))**2.d0
!     &            *(DX(1)*beta*DYF(J)*2.d0*DZ(1)*beta)**(2.d0/3.d0)

      END DO
! Extend to ghost cells 
      DELTA_YF(0)=DELTA_YF(1)
      DELTA_YF(NY+1)=DELTA_YF(NY)      

      DO J=1,NY+1
! At GY points:
! Constant Smagorinsky
        DELTA_Y(J)=-C_SMAG**2.d0
     &        *(DX(1)*beta*DY(J)*2.d0*DZ(1)*beta)**(2.d0/3.d0)
! Wall Damping
!        DELTA_Y(J)=
!     &    -1.d0*(0.1d0*(1.d0-exp((-GY(J)/(NU*25.d0))**3.d0)))**2.d0
!     &            *(DX(1)*beta*DY(J)*2.d0*DZ(1)*beta)**(2.d0/3.d0)
      END DO

! Get the eddy diffusivity at GY points
! NU_T = (C_S^2 * DELTA^2)*|S| With |S| interpolated to GY points
      DO J=J1,J2
        DO K=0,NZM
          DO I=0,NXM
            KAPPA_T(I,K,J,N)=-1.d0*DELTA_Y(J)
     &       *(S1(I,K,J)*DYF(J-1)+S1(I,K,J-1)*DYF(J))/(2.d0*DY(J))
          END DO
        END DO
      END DO

! Now, compute TAU, store in the corresponging Sij
      DO K=0,TNKZ
        DO I=0,NKX
          DO J=1,NY
            CSij(I,K,J,1)=DELTA_YF(J)*CSij(I,K,J,1)
            CSij(I,K,J,3)=DELTA_YF(J)*CSij(I,K,J,3)
          END DO
        END DO
      END DO

! tau_ij is now contained in CSij in Fourier space

! Convert the scalar to physical space
      call FFT_XZ_TO_PHYSICAL(CTH(0,0,0,N),TH(0,0,0,N),0,NY+1)

      else if ((LES_MODEL_TYPE.EQ.2).or.(LES_MODEL_TYPE.eq.3)) then
! Here, use a dynamic smagorinsky model
! Note, there is no scale similar model for the scalar,
! so model type choice 2 and 3 are identical for the scalar equation

! Compute the filter width
      DO J=0,NY+1
! At GYF points:
        DELTA_YF(J)=(beta*DX(1)*DYF(J)*beta*DZ(1))**(1.d0/3.d0)
!        DELTA_YF(J)=sqrt((beta*DX(1))**2.d0+(DYF(J)*2.d0)**2.d0
!     &      +(beta*DZ(1))**2.d0)
      END DO
! Extend to ghost cells
!      DELTA_YF(0)=DELTA_YF(1)
!      DELTA_YF(NY+1)=DELTA_YF(NY)

      DO J=1,NY+1
! At GY points:
       DELTA_Y(J)=beta*(beta*DX(1)*DY(J)*beta*DZ(1))**(1.d0/3.d0)
!       DELTA_Y(J)=sqrt((beta*DX(1))**2.d0+(DY(J)*2.d0)**2.d0
!     &          +(beta*DZ(1))**2.d0)
      END DO

! We need to calculate the components of C, the dynamic coefficient
! C_DYN_TH will be defined at GYF points

! Compute the scalar gradient, store in Sij(:,:,:,1..3)
      call compute_scalar_grad(n)

! Convert the scalar to physical space
      CALL FFT_XZ_TO_PHYSICAL(CTH(0,0,-1,N),TH(0,0,-1,N),0,NY+1)
      CALL FFT_XZ_TO_PHYSICAL(CTH(0,0,NY+1,N),TH(0,0,NY+1,N),0,1)

! Compute C_DYN_TH only every x # of timesteps
      if ((MOD(TIME_STEP,15).eq.0).OR.FIRST_TIME) THEN

! Store TH in Sij(:,:,:,4) and apply the test filter
      do j=0,NY+1
        do k=0,NZM
          do i=0,NXM
            Sij(i,k,j,4)=TH(i,k,j,n)
          end do
        end do
      end do
      call les_filter_chan(Sij(0,0,0,4),0,NY+1,2)

! Zero C_DYN_TH
      do j=0,NY+1
        C_DYN_TH(j,N)=0.d0
        denominator_sum(j)=0.d0 
        numerator_sum(j)  =0.d0
      end do

! Do over all non-repeating components of the scalar gradient
      do ij=1,3

! Zero the numerator and denominator:
        do j=0,NY+1
          do k=0,NZ+1
            do i=0,NX+1
              numerator(i,k,j)=0.d0
              denominator(i,k,j)=0.d0
            end do
          end do
        end do
        
! First, compute Mij

       do k=0,NZM
         do i=0,NXM
           do j=0,min(NY+1,J2+1)
             if (ij.eq.2) then
! Sij is defined at GY points, interpolate to GYF points           
               temp(i,k,j)=0.5d0*(Sij(i,k,j+1,ij)+Sij(i,k,j,ij))
             else
! Sij is defined at GYF points, no interpolation needed
               temp(i,k,j)=Sij(i,k,j,ij)
             end if
! Define temp at ghost points since it will be filtered
!             temp(i,k,0)=temp(i,k,1)
!             temp(i,k,NY+1)=temp(i,k,NY)
           end do
         end do
       end do
! Filter temp
       call les_filter_chan(temp,0,min(NY+1,J2+1),2)
! Multiply by |S| filtered        
       do j=0,min(NY+1,J2+1)
         do k=0,NZM
           do i=0,NXM
             temp(i,k,j)=temp(i,k,j)*(alpha*DELTA_YF(j))**2.d0
     &                                    *S_2BAR(i,k,j) 
           end do
         end do
       end do
! Get second term of Mij
       if (ij.eq.2) then
! Sij is defined at GY points, interpolate to GYF points
         do k=0,NZM
           do i=0,NXM
             do j=0,min(NY+1,J2+1)
               Mij(i,k,j)=DELTA_YF(j)**2.d0*S1(i,k,j)
     &               *0.5d0*(Sij(i,k,j+1,ij)+Sij(i,k,j,ij))
             end do
! Define Mij at ghost cells since it will be filtered
!             Mij(i,k,0)=Mij(i,k,1)
!             Mij(i,k,NY+1)=Mij(i,k,NY)
           end do
         end do
       else
! Sij is defined at GYF points, no interpolation needed
         do i=0,NXM
           do k=0,NZM 
             do j=0,min(NY+1,J2+1)
               Mij(i,k,j)=DELTA_YF(j)**2.d0*S1(i,k,j)*Sij(i,k,j,ij)
             end do
! Define Mij at ghost cells
!             Mij(i,k,0)=Mij(i,k,1)
!             Mij(i,k,NY+1)=Mij(i,k,NY)
           end do
         end do
       end if

! Filter Mij
       call les_filter_chan(Mij,0,min(NY+1,J2+1),2)
 
! Add the second term of Mij stored in temp
       Mij=temp-Mij    
! Now, compute Lij and add Lij*Mij to the numerator
! temp=Ui*Uj:
       SELECT CASE (ij)
       CASE(1)
         do j=0,min(NY+1,J2+1)
           do k=0,NZM
             do i=0,NXM
               temp(i,k,j)=U1(i,k,j)*TH(i,k,j,n)
             end do
           end do 
         end do
       CASE(2)
         do j=0,min(NY+1,J2+1)
           do k=0,NZM
             do i=0,NXM
               temp(i,k,j)=TH(i,k,j,n)*0.5d0*(U2(i,k,j+1)+U2(i,k,j))
             end do
           end do
         end do
! Get estimate at ghost cells, Use U2 at GY ghost cell directly
!         do k=0,NZM
!           do i=0,NXM
!             temp(i,k,0)=TH(i,k,0,n)*U2(i,k,1)
!             temp(i,k,NY+1)=TH(i,k,NY+1,n)*U2(i,k,NY+1)
!           end do 
!         end do
       CASE(3)
         do j=0,min(NY+1,J2+1)
           do k=0,NZM
             do i=0,NXM
               temp(i,k,j)=U3(i,k,j)*TH(i,k,j,n)
             end do
           end do 
         end do
       END SELECT
! Filter temp
       call les_filter_chan(temp,0,NY+1,2)
! Add Lij*Mij to numerator
! Recall that Sij(:,:,:,4) holds TH_2BAR
       do j=0,min(NY+1,J2+1)
         do k=0,NZM
           do i=0,NXM 
             numerator(i,k,j)=Mij(i,k,j)
     &              *(temp(i,k,j)-Sij(i,k,j,4)*U_BAR_TIL(i,k,j,ij))
           end do
         end do
       end do

! Now, the denominator for this ij  
       do j=0,NY+1
         do k=0,NZM
           do i=0,NXM
             denominator(i,k,j)=Mij(i,k,j)*Mij(i,k,j)
           end do
         end do
       end do

! Get plane average of numerator and denominator, add to C
! Note, since both the numerator and denominator are averaged, there
! is not need to divide the sum by NX*NZ

      do j=0,NY+1
        denominator_sum(j)=denominator_sum(j)+
     &    SUM(denominator(0:NXM,0:NZM,j))
        numerator_sum(j)=numerator_sum(j)+SUM(numerator(0:NXM,0:NZM,j))
      end do

! End to ij
      end do


       do j=0,NY+1
        if (denominator_sum(j) .ne. 0.) then
          C_DYN_TH(j,n)=-0.5d0*numerator_sum(j)/denominator_sum(j)
        else
          C_DYN_TH(j,n)=0.d0
        endif
       enddo

! We are now done with the dynamic procedure to calculate C

! If C_DYN_TH < 0 at any level, set C_DYN_TH=0 for numerical stability
      do j=0,NY+1
        if (C_DYN_TH(j,n).lt.0) C_DYN_TH(j,n)=0.d0
      end do

! End if compute C_DYN_TH
       END IF

! Get the eddy diffusivity at GY points
! KAPPA_T = C_DYN_TH * (DELTA^2)*|S|
! Use exact second order interpolation to get C_DYN_TH and S1 interpolated to
! GY points
c      DO J=J1,J2
      DO J=J1i,min(NY+1,J2+1)
        DO K=0,NZM
          DO I=0,NXM
            KAPPA_T(I,K,J,N)=
     &        (DYF(j-1)*C_DYN_TH(j,n)+DYF(j)*C_DYN_TH(j-1,n))
     &            /(2.d0*DY(j))
     &        *DELTA_Y(j)**2.d0
     &        *(DYF(j-1)*S1(i,k,j)+DYF(j)*S1(i,k,j-1))/(2.d0*DY(j))
          END DO
        END DO
      END DO

c      write(400,*) time_step, sum(KAPPA_T(0:NXM,0:NZM,1,1)), 10000
c      write(410,*) time_step, sum(KAPPA_T(0:NXM,0:NZM,NY,1)), 20000

c      DO J=J2+2,NY+1
c        DO K=0,NZM
c          DO I=0,NXM
C            KAPPA_T(I,K,J,N)=0.d0 
c          ENDDO
c        ENDDO
c      ENDDO 

! At this point we have C_DYN_TH and dTH/dx_i (stored in Sij(:,:,:,1...3)
! Calculate lambda_i in physical space, stored in Sij(:,:,:,1..3)

      do ij=1,3
! Dynamic Smagorinsky model, no scale similar part
! Here, Sij is defined at GYF points 
      do k=0,NZM
        do i=0,NXM
          do j=0,min(NY+1,J2+1)
            Sij(i,k,j,ij)=
     &        -C_DYN_TH(j,n)*DELTA_YF(j)**2.d0*S1(i,k,j)*Sij(i,k,j,ij)
          end do 
! Get TAUij at ghost cells for computing derivative
C          Sij(i,k,0,ij)=Sij(i,k,1,ij)
C          Sij(i,k,NY+1,ij)=Sij(i,k,NY,ij)
        end do
       end do

! Convert TAUij, now stored in Sij to Fourier space
       call FFT_XZ_TO_FOURIER(Sij(0,0,0,ij),CSij(0,0,0,ij),0
     &          ,min(NY+1,J2+1))

! End do ij
       end do  

! End if LES_MODEL_TYPE dynamic Smagorinsky or Dynamic model
      else
        stop
        write(6,*) 'Error, unsupported LES_MODEL_TYPE chosen'
      end if

! Now, add the subgrid scale forcing to CFi
! (This includes the subgrid scale stress as an explicit R-K term
! Include only CSij terms 1 and 3 since term 2 is accounted for
! as an implicit eddy diffusivity through KAPPA_T

!      DO J=J2+1,NY+1
!        DO K=0,TNKZ
!          DO I=0,NKX
!            DO ij =1,6
!             CSij(I,K,J,ij) = 0.d0
!            ENDDO
!           ENDDO
!         ENDDO
!       ENDDO

      DO K=0,TNKZ
        DO I=0,NKX
          DO J=J1,J2
            CFTH(I,K,J,n)=CFTH(I,K,J,n)
     &                -CIKX(I)*CSij(I,K,J,1)
     &                -CIKZ(K)*CSij(I,K,J,3)
          END DO
       END DO
      END DO

! Periodically, output mean quantities
      IF ((MOD(TIME_STEP,SAVE_STATS_INT).EQ.0).AND.(RK_STEP.EQ.1)) THEN
c       open(43,file='mean_les_th.txt',form='formatted',status='unknown')
       open(43,file=FNAME_TH,form='formatted',status='unknown',
     & position='append')
       write(43,*) TIME_STEP,TIME,DELTA_T
       do j=1,NY
          write(43,420) j,GYF(J)
     &    ,(dble(CSij(0,0,J,ij)),ij=1,3),
     &    SUM(KAPPA_T(0:NXM,0:NZM,J,N))/DBLE(NX*NZ),
     &       C_DYN_TH(J,n)*DELTA_YF(J)**2.d0
       end do
      close(43)
      END IF
420     format(I3,' ',6(F20.9,' '))

      RETURN
      END


      subroutine compute_scalar_grad(n)
C This subroutine computes dTH/dx_i for the filtered scalar field
C The input velocity field should be in fourier space in the periodic
C directions.
C For use in the LES model in channel flow (2 periodic directions)
C Store in Sij(:,:,:,1..3)
      include 'header_les'

      integer I,J,K,ij,n

      DO K=0,TNKZ
        DO I=0,NKX
          DO J=0,NY+1
            CSij(I,K,J,1)=CIKX(I)*CTH(I,K,J,n)
            CSij(I,K,J,3)=CIKZ(K)*CTH(I,K,J,n)
          END DO
          DO J=0,NY+2
! Define the vertical gradient at GY points
            CSij(I,K,J,2)=(CTH(I,K,J,n)-CTH(I,K,J-1,n))/DY(j)
          END DO
        END DO
      END DO

! Convert the scalar gradients to physical space
      do ij=1,3
        call FFT_XZ_TO_PHYSICAL(CSij(0,0,0,ij),Sij(0,0,0,ij),0,NY)
        call FFT_XZ_TO_PHYSICAL(CSij(0,0,NY+1,ij),Sij(0,0,NY+1,ij),0,1)
      end do

! We now have dTH/dx_i in Physical space
      RETURN
      END







