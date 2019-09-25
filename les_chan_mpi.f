      subroutine les_chan
C This subroutine models the terms owing to the subgrid scale stress
C if the computation is to be treated as an LES not a DNS
C This subroutine should be called when the velocity is in fourier space 
C in the periodic directions, on output, the velocity will be 
C in physical space.
C It is assumed that the test filter and the LES filter are performed
C by the same operation
C On output S1 should contain |S| which may be used again in les_chan_th
C if for the subgrid scalar dissipation

      include 'header_les'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'

      integer i,j,k,l,m,ij
      CHARACTER*35 FNAME, FNAME_tke
      real*8 S1_mean(1:NY)
      real*8 NU_T_mean(1:NY)
      real*8 tke_sgs_t(0:NY+1),tke_sgs_p(0:NY+1),
     &       tke_sgs_diss(0:NY+1)

! Variables for Dynamic Smagorinsky model:
      real*8 C_SMAG
      parameter (C_SMAG=0.13d0)
      real*8 DELTA_Y(0:NY+1) 
      real*8 alpha,beta 
      real*8 denominator_sum(0:NY+1),
     &       numerator_sum(0:NY+1)  
! Array to store the velocity index for each component of the strain rate tensor
      integer U_index1(6)
      integer U_index2(6)

! Here, alpha is the test/LES filter width ratio
      parameter (alpha= 1.6 ) !2.44950)
! beta is the LES/grid filter width ratio
      parameter (beta= 1.d0)

! Set the velocity index for each component of the stress tensor
      U_index1(1)=1
      U_index2(1)=1

      U_index1(2)=2
      U_index2(2)=2 

      U_index1(3)=3
      U_index2(3)=3

      U_index1(4)=1
      U_index2(4)=2

      U_index1(5)=1
      U_index2(5)=3

      U_index1(6)=2
      U_index2(6)=3

       I = RANK+1

           FNAME ='results/mean_les_'
     &        //CHAR(MOD(I,100)/10+48)
     &        //CHAR(MOD(I,10)+48) // '.txt'

           FNAME_TKE ='results/tke_les_'
     &        //CHAR(MOD(I,100)/10+48)
     &        //CHAR(MOD(I,10)+48) // '.txt'

! When using a Near-wall model, don't use LES at the wall
C      IF ((U_BC_YMIN.EQ.3).or.(W_BC_YMIN.EQ.3)) then
C         J1=2
C      ELSE
C         J1=JSTART
C      END IF

C      IF ((U_BC_YMAX.EQ.3).or.(W_BC_YMAX.EQ.3)) then
C         J2=NY-1
C      ELSE
C         J2=JEND
C      END IF

c       write(400,*) time_step, CTH(0,0,1,1), 10000,RK_STEP
c       write(410,*) time_step, CTH(0,0,NY,1), 20000,RK_STEP

! First, for all models, apply boundary conditions to the velocity field
! (fill ghost cells) to ensure accurate calculation of gradients
C Apply Boundary conditions to velocity field
      IF (USE_MPI) THEN
        J1  = JSTART
        J2  = JEND
        J1i= 1
        J2e=NY+1
        CALL APPLY_BC_VEL_MPI
      ELSE
        J1i=2
        J2e=NY
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

c      write(333,*) J1i, J2e

      if (LES_MODEL_TYPE.EQ.1) then
!     Constant Smagorinsky model

! First, compute the rate of strain tensor S_ij

      call compute_strain_chan

! Compute |S| at GYF points, store in S1
! Interpolation to GYF points is easy since by definition
! GYF points are exactly midway between neighboring GY points
      DO J=0,NY+1
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,J)=SQRT(
     &                2.d0*Sij(I,K,J,1)**2.d0
     &               +4.d0*(0.5d0*(Sij(I,K,J+1,4)+Sij(I,K,J,4)))**2.d0
     &               +4.d0*Sij(I,K,J,5)**2.d0
     &               +2.d0*Sij(I,K,J,2)**2.d0
     &               +4.d0*(0.5d0*(Sij(I,K,J+1,6)+Sij(I,K,J,4)))**2.d0 
     &               +2.d0*Sij(I,K,J,3)**2.d0 )
          END DO
        END DO
      END DO
! Extend |S| to ghost cells
      DO K=0,NZM
        DO I=0,NXM
          S1(I,K,0)=S1(I,K,1)
          S1(I,K,NY+1)=S1(I,K,NY)
        END DO
      END DO


! Now, compute |S|*S_ij, storing in Sij
! First compute at GYF points 
      DO J=1,NY
        DO K=0,NZM
          DO I=0,NXM
            Sij(I,K,J,1)=S1(I,K,J)*Sij(I,K,J,1)
            Sij(I,K,J,5)=S1(I,K,J)*Sij(I,K,J,5)
! CSij(:,:,:,2) is added through an implicit eddy viscosity
            Sij(I,K,J,2)=0.d0
            Sij(I,K,J,3)=S1(I,K,J)*Sij(I,K,J,3)
          END DO
        END DO
      END DO
! Now, compute at GY points, interpolating |S|
      DO J=1,NY+1
        DO K=0,NZM
          DO I=0,NXM
! |S| interpolated to GY point 
            TEMP(I,K,J)=(S1(I,K,J)*DYF(j-1)+S1(I,K,J-1)*DYF(j))
     &                  /(2.d0*DY(j))
! The terms dU1/dy and dU3/dy in CSij(:,:,:,4) and CSij(:,:,:,6) respectively
! are subtracted out from Sij here since they are treated implicitly
! in eddy viscosity terms
            Sij(I,K,J,4)=TEMP(I,K,J)
     &        *(Sij(I,K,J,4)-0.5*(CU1(I,K,J)-CU1(I,K,J-1))/DY(j))
            Sij(I,K,J,6)=TEMP(I,K,J)
     &        *(Sij(I,K,J,6)-0.5*(CU3(I,K,J)-CU3(I,K,J-1))/DY(j))
          END DO
        END DO
      END DO



! We now have |S|*S_ij stored in Sij in Physical space

! Convert |S|*S_ij to Fourier space
      do ij=1,6 
        CALL FFT_XZ_TO_FOURIER(Sij(0,0,0,ij),CSij(0,0,0,ij),0,NY+1)
      end do
! Compute the filter lengthscale
! Absorb -2.d0*C_SMAG**2.d0 here for effienciency
      DO J=1,NY
! At GYF points:
! Constant Smagorinsky
!        DELTA_YF(J)=-2.d0*C_SMAG**2.d0
!     &     *(DX(1)*beta*DYF(J)*2.d0*DZ(1)*beta)**(2.d0/3.d0)
! Wall Damping
        DELTA_YF(J)=
     &    -2.d0*(0.1d0*(1.d0-exp((-GYF(J)/(NU*25.d0))**3.d0)))**2.d0
     &            *(DX(1)*beta*DYF(J)*2.d0*DZ(1)*beta)**(2.d0/3.d0)

      END DO
! Extend to ghost cells 
      DELTA_YF(0)=DELTA_YF(1)
      DELTA_YF(NY+1)=DELTA_YF(NY)      

      DO J=1,NY+1
! At GY points:
! Constant Smagorinsky
!        DELTA_Y(J)=-2.d0*C_SMAG**2.d0
!     &        *(DX(1)*beta*DY(J)*2.d0*DZ(1)*beta)**(2.d0/3.d0)
! Wall Damping
        DELTA_Y(J)=
     &    -2.d0*(0.1d0*(1.d0-exp((-GY(J)/(NU*25.d0))**3.d0)))**2.d0
     &            *(DX(1)*beta*DY(J)*2.d0*DZ(1)*beta)**(2.d0/3.d0)
      END DO

! Get the eddy viscosity at GY points
! NU_T = (C_S^2 * DELTA^2)*|S|

      DO J=1,NY+1
        DO K=0,NZM
          DO I=0,NXM
            NU_T(I,K,J)=-0.5d0*DELTA_Y(J)*TEMP(I,K,J)
          END DO
        END DO
      END DO

! Now, compute TAU, store in the corresponging Sij
      DO J=1,NY
        DO K=0,TNKZ
          DO I=0,NKX
            CSij(I,K,J,1)=DELTA_YF(J)*CSij(I,K,J,1)
            CSij(I,K,J,5)=DELTA_YF(J)*CSij(I,K,J,5)
! CSij(:,:,:,2) is added through an implicit eddy viscosity
!            CSij(I,K,J,2)=DELTA_YF(J)*CSij(I,K,J,2)
            CSij(I,K,J,3)=DELTA_YF(J)*CSij(I,K,J,3)
          END DO
        END DO
      END DO
      DO J=1,NY+1 
        DO K=0,TNKZ
          DO I=0,NKX
            CSij(I,K,J,4)=DELTA_Y(J)*CSij(I,K,J,4)
            CSij(I,K,J,6)=DELTA_Y(J)*CSij(I,K,J,6)
          END DO
        END DO
      END DO

! tau_ij is now contained in CSij in Fourier space

! Convert the velocity to physical space
      call FFT_XZ_TO_PHYSICAL(CU1(0,0,0),U1(0,0,0),0,NY+1)
      call FFT_XZ_TO_PHYSICAL(CU2(0,0,0),U2(0,0,0),0,NY+1)
      call FFT_XZ_TO_PHYSICAL(CU3(0,0,0),U3(0,0,0),0,NY+1)

      else if ((LES_MODEL_TYPE.EQ.2).or.(LES_MODEL_TYPE.eq.3)) then
! Here, use a dynamic smagorinsky model with or without scale similar part

! Compute the filter width
      DO J=0,NY+1
! At GYF points:
        DELTA_YF(J)=(beta*DX(1)*DYF(J)*beta*DZ(1))**(1.d0/3.d0)
!        DELTA_YF(J)=sqrt((beta*DX(1))**2.d0+(DYF(J)*2.d0)**2.d0
!     &      +(beta*DZ(1))**2.d0)
      END DO

      IF (USE_MPI) THEN
          CALL  GHOST_ZERO_CHAN_MPI
      ENDIF

! Extend to ghost cells
!      DELTA_YF(0)=DELTA_YF(1)
!      DELTA_YF(NY+1)=DELTA_YF(NY)

      DO J=0,NY+1
! At GY points:
        DELTA_Y(J)=(beta*DX(1)*DY(J)*beta*DZ(1))**(1.d0/3.d0)
!       DELTA_Y(J)=sqrt((beta*DX(1))**2.d0+(DY(J)*2.d0)**2.d0
!     &          +(beta*DZ(1))**2.d0)
      END DO

! We need to calculate the components of C, the dynamic coefficient
! C_DYN will be defined at GYF points

! Compute the rate of strain tensor, store in Sij
      call compute_strain_chan

 
! Compute |S| at GYF points, store in S1
! Interpolation to GYF points is easy since by definition
! GYF points are exactly midway between neighboring GY points
      DO J=0,NY+1
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,J)=SQRT(
     &                2.d0*Sij(I,K,J,1)**2.d0
     &               +4.d0*(0.5d0*(Sij(I,K,J+1,4)+Sij(I,K,J,4)))**2.d0
     &               +4.d0*Sij(I,K,J,5)**2.d0
     &               +2.d0*Sij(I,K,J,2)**2.d0
     &               +4.d0*(0.5d0*(Sij(I,K,J+1,6)+Sij(I,K,J,6)))**2.d0 
     &               +2.d0*Sij(I,K,J,3)**2.d0 )
          END DO
        END DO
      END DO

! Extend |S| to ghost cells
!      DO K=0,NZM
!        DO I=0,NXM
!          S1(I,K,0)=S1(I,K,1)
!          S1(I,K,NY+1)=S1(I,K,NY)
!        END DO
!      END DO
        

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!      Storing for tke les
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      DO ij=1,6
       DO J=0,NY+1
         DO K=0,NZM
          DO I=0,NXM
            St_rij(I,K,J,ij)=Sij(I,K,J,ij)
          ENDDO
         ENDDO
       ENDDO
      ENDDO


      DO J=0,NY+1
       U1_BAR(j) = CU1(0,0,j)
       U2_BAR(j) = CU2(0,0,j)
       U3_BAR(j) = CU3(0,0,j)
      ENDDO
! Convert Ui to physical space
      CALL FFT_XZ_TO_PHYSICAL(CU1(0,0,-1),U1(0,0,-1),0,NY+1)
      CALL FFT_XZ_TO_PHYSICAL(CU2(0,0,-1),U2(0,0,-1),0,NY+1)
      CALL FFT_XZ_TO_PHYSICAL(CU3(0,0,-1),U3(0,0,-1),0,NY+1)

      CALL FFT_XZ_TO_PHYSICAL(CU1(0,0,NY+1),U1(0,0,NY+1),0,1)
      CALL FFT_XZ_TO_PHYSICAL(CU2(0,0,NY+1),U2(0,0,NY+1),0,1)
      CALL FFT_XZ_TO_PHYSICAL(CU3(0,0,NY+1),U3(0,0,NY+1),0,1)       

      


! Apply the filter to the LES velocity and save
      do j=0,NY+1
        do k=0,NZM
          do i=0,NXM
            U_BAR_TIL(i,k,j,1)=U1(i,k,j)
            U_BAR_TIL(i,k,j,3)=U3(i,k,j)
          end do
        end do
      end do
      do j=0,NY+1
        do k=0,NZM
          do i=0,NXM
! Interpolate U2 to GYF points
            U_BAR_TIL(i,k,j,2)=0.5d0*(U2(i,k,j+1)+U2(i,k,j))
          end do
        end do
      end do

!      do k=0,NZM
!        do i=0,NXM
! Extrapolate to GYF ghost cells
!       U_BAR_TIL(i,k,0,2)=U2(i,k,1)*(1.d0-DY(1)/(2.d0*DYF(1)))
!     &                     -U2(i,k,2)*DY(1)/(2.d0*DYF(1))
!       U_BAR_TIL(i,k,NY+1,2)=U2(i,k,NY+1)*(1.d0+DY(NY+1)/(2.d0*DYF(NY)))
!     &                     -U2(i,k,NY)*DY(NY+1)/(2.d0*DYF(NY)) 
!        end do
!      end do       

       

      
      if (les_model_type.eq.3) then
! storing in the U_2BAR
       do j=0,NY+1
         do k=0,NZM
           do i=0,NXM
             do ij =1,3
               U_2BAR(i,k,j,ij) = U_BAR_TIL(i,k,j,ij)
             end do
           end do
         end do
       end do

       do i=1,3
! application of grid filter i.e. filter_type = 1
        call les_filter_chan(U_2BAR(0,0,0,i),0,NY+1,1)
       enddo
      endif
   
! Now, filter the velocity
      do i=1,3
! application of test filter i.e. filter_type = 2
        call les_filter_chan(U_BAR_TIL(0,0,0,i),0,NY+1,2)
      end do

! Compute C_DYN only every x # of timesteps
      if (((MOD(TIME_STEP,15).eq.0).AND.(RK_STEP.eq.1))
     &       .or.FIRST_TIME) THEN

! Filter |S| and store in S_2BAR
      DO J=0,NY+1 
        DO K=0,NZM
          DO I=0,NXM
            S_2BAR(I,K,J)=S1(I,K,J)
          END DO
        END DO
      END DO
! Test filtering operation filter type = 2
      call les_filter_chan(S_2BAR,0,NY+1,2)

! Save a copy of the velocity which will be filtered twice more 

      if (les_model_type.eq.3) then 
! Do only if a scale similar part is needed
        do ij=1,3  
        do j=0,NY+1
          do k=0,NZ+1
            do i=0,NX+1
              U_4BAR(i,k,j,ij)=U_BAR_TIL(i,k,j,ij)
            end do
          end do
        end do
        end do

        do i=1,3
! application of test and bar filter togather i.e. two filtering
! operation first filter type = 1, second filter type = 2
          call les_filter_chan(U_4BAR(0,0,0,i),0,NY+1,1)
          call les_filter_chan(U_4BAR(0,0,0,i),0,NY+1,2)
        end do
      end if


! Zero C_DYN
      do j=0,NY+1
        C_DYN(j)=0.d0
      end do

! The prep. work is now done, we are ready to start the algorithm to compute
! the dynamic model coefficient

      DO j =0,NY+1
       denominator_sum(j) = 0.d0
       numerator_sum(j)   = 0.d0
      ENDDO
 
! Do over all non-repeating components of the stress tensor
      do ij=1,6

! Here        ij=1 -> l=1,m=1
!             ij=2 -> l=2,m=2
!             ij=3 -> l=3,m=3
!             ij=4 -> l=1,m=2
!             ij=5 -> l=1,m=3
!             ij=6 -> l=2,m=3       
  
! Zero the numerator and denominator:
        do j=0,NY+1
          do k=0,NZ+1
            do i=0,NX+1
              numerator(i,k,j)=0.d0
              denominator(i,k,j)=0.d0
            end do
          end do
        end do

! cross is used to multiply by two to include the contribution from the
! other symmetric term if we are dealing with a cross-term
        if (ij.le.3) then 
! We are computing a diagonal term
          cross=1.d0
        else
! We are computing a cross term
          cross=2.d0 
        end if 

! First, compute Mij

       if ((ij.eq.4).or.(ij.eq.6)) then
       do j=0,NY+1 
         do k=0,NZM
           do i=0,NXM
! Sij is defined at GY points, interpolate to GYF points           
               temp(i,k,j)=0.5d0*(Sij(i,k,j+1,ij)+Sij(i,k,j,ij))
           end do
         end do
       end do
       else
       do j=0,NY+1 
         do k=0,NZM
           do i=0,NXM
! Sij is defined at GYF points, no interpolation needed
             temp(i,k,j)=Sij(i,k,j,ij)
           end do
         end do
        end do
        end if
! Define temp at ghost points since it will be filtered
!        do k=0,NZM
!          do i=0,NXM
!             temp(i,k,0)=temp(i,k,1)
!             temp(i,k,NY+1)=temp(i,k,NY)
!           end do
!         end do
! Filter temp test filter operation filter type = 2
       call les_filter_chan(temp,0,NY+1,2)
! Multiply by |S| filtered        
       do j=0,NY+1
         do k=0,NZM
           do i=0,NXM
             temp(i,k,j)=temp(i,k,j)*(alpha*DELTA_YF(j))**2.d0
     &                                    *S_2BAR(i,k,j) 
           end do
         end do
       end do

      
 
! Get second term of Mij
       if ((ij.eq.4).or.(ij.eq.6)) then
! Sij is defined at GY points, interpolate to GYF points
         do j=0,NY+1
           do k=0,NZM
             do i=0,NXM
               Mij(i,k,j)=DELTA_YF(j)**2.d0*S1(i,k,j)
     &               *0.5d0*(Sij(i,k,j+1,ij)+Sij(i,k,j,ij))
             end do
           end do
         end do
! Define Mij at ghost cells since it will be filtered
!         do k=0,NZM
!           do i=0,NXM
!             Mij(i,k,0)=Mij(i,k,1)
!             Mij(i,k,NY+1)=Mij(i,k,NY)
!           end do
!         end do
       else
! Sij is defined at GYF points, no interpolation needed
         do j=0,NY+1
           do i=0,NXM
             do k=0,NZM 
               Mij(i,k,j)=DELTA_YF(j)**2.d0*S1(i,k,j)*Sij(i,k,j,ij)
             end do
           end do
         end do
! Define Mij at ghost cells
!         do i=0,NXM
!           do k=0,NZM 
!             Mij(i,k,0)=Mij(i,k,1)
!             Mij(i,k,NY+1)=Mij(i,k,NY)
!           end do
!         end do
       end if

! Filter Mij test filter operation filter type = 2
       call les_filter_chan(Mij,0,NY+1,2)
 
! Add the second term of Mij stored in temp
       do j=0,NY+1
         do k=0,NZM
           do i=0,NXM
             Mij(i,k,j)=temp(i,k,j)-Mij(i,k,j)    
           end do
         end do
       end do
! Now, compute Lij and add Lij*Mij to the numerator
! temp=Ui*Uj:
       SELECT CASE (ij)
       CASE(1)
         do j=0,NY+1
           do k=0,NZM
             do i=0,NXM
               temp(i,k,j)=U1(i,k,j)*U1(i,k,j)
             end do
           end do 
         end do
       CASE(2)
         do j=0,NY+1
           do k=0,NZM
             do i=0,NXM
               temp(i,k,j)=0.25d0*(U2(i,k,j+1)+U2(i,k,j))**2.d0
             end do
           end do
         end do
! Get estimate at ghost cells, Use U2 at GY ghost cell directly
!         do k=0,NZM
!           do i=0,NXM
!               temp(i,k,0)=U2(i,k,1)**2.d0
!               temp(i,k,NY+1)=U2(i,k,NY+1)**2.d0
!           end do 
!         end do
       CASE(3)
         do j=0,NY+1
           do k=0,NZM
             do i=0,NXM
               temp(i,k,j)=U3(i,k,j)*U3(i,k,j)
             end do
           end do 
         end do
       CASE(4) 
         do j=0,NY+1
           do k=0,NZM
             do i=0,NXM
               temp(i,k,j)=U1(i,k,j)*0.5d0*(U2(i,k,j+1)+U2(i,k,j))
             end do
           end do
         end do
! Get estimate at ghost cells, Use U2 at GY ghost cell directly
!         do k=0,NZM
!           do i=0,NXM
!             temp(i,k,0)=U1(i,k,0)*U2(i,k,1)
!             temp(i,k,NY+1)=U1(i,k,NY+1)*U2(i,k,NY+1)
!           end do 
!         end do
       CASE(5)
         do j=0,NY+1
           do k=0,NZM
             do i=0,NXM
               temp(i,k,j)=U1(i,k,j)*U3(i,k,j)
             end do
           end do 
         end do
       CASE(6) 
         do j=0,NY+1
           do k=0,NZM
             do i=0,NXM
               temp(i,k,j)=0.5d0*(U2(i,k,j+1)+U2(i,k,j))*U3(i,k,j)
             end do
           end do
         end do
! Get estimate at ghost cells, Use U2 at GY ghost cell directly
!         do k=0,NZM
!           do i=0,NXM
!             temp(i,k,0)=U3(i,k,0)*U2(i,k,1)
!             temp(i,k,NY+1)=U3(i,k,NY+1)*U2(i,k,NY+1)
!           end do 
!         end do
       END SELECT

       do k=0,NZM
           do i=0,NXM
             do j=0,NY+1
              temp_1(i,k,j) = temp(i,k,j)
             end do
           end do
       end do
! Filter temp: test filter operation i.e. filter type = 2
       call les_filter_chan(temp,0,NY+1,2)
! Add Lij*Mij to numerator
       do j=0,NY+1
         do k=0,NZM
           do i=0,NXM 
             numerator(i,k,j)=Mij(i,k,j)
     &        *(temp(i,k,j)
     &     -U_BAR_TIL(i,k,j,U_index1(ij))*U_BAR_TIL(i,k,j,U_index2(ij)))
           end do
         end do
       end do

       if (LES_MODEL_TYPE.eq.3) then
! If mixed model, include the scale similar part  
! Add Nij*Mij to the numerator piece-by-piece
! Term3 

! similarly grid filter operation is done on temp i.e. filter type = 1
         call les_filter_chan(temp_1,0,NY+1,1)
! Filter temp_1 test filter operation filter type = 2
         call les_filter_chan(temp_1,0,NY+1,2)
         do j=0,NY+1
           do k=0,NZM
             do i=0,NXM
              numerator(i,k,j)=numerator(i,k,j)+temp_1(i,k,j)*Mij(i,k,j)
             end do
           end do
         end do
! Term 4
         do j=0,NY+1
           do k=0,NZM
             do i=0,NXM
               temp(i,k,j)=U_2BAR(i,k,j,U_index1(ij))
     &                    *U_2BAR(i,k,j,U_index2(ij))
               temp_1(i,k,j) = U_BAR_TIL(i,k,j,U_index1(ij))
     &                    *U_BAR_TIL(i,k,j,U_index2(ij))
             end do
           end do
         end do
! Filter temp test filter operation filter type = 2
         call les_filter_chan(temp,0,NY+1,2)
         do j=0,NY+1
           do k=0,NZM
             do i=0,NXM
               numerator(i,k,j)=numerator(i,k,j)-temp(i,k,j)*Mij(i,k,j)
             end do
           end do
         end do 
! Terms 1 and 2
! Filter temp_1 grid filter operation filter type = 1
         call les_filter_chan(temp_1,0,NY+1,1)
! Filter temp test filter operation filter type = 2
         call les_filter_chan(temp_1,0,NY+1,2)
         do j=0,NY+1
           do k=0,NZM
             do i=0,NXM
              numerator(i,k,j)=numerator(i,k,j)
     &        -(temp_1(i,k,j)
     &         -U_4BAR(i,k,j,U_index1(ij))*U_4BAR(i,k,j,U_index2(ij)))
     &           *Mij(i,k,j)
             end do
           end do
         end do
       end if
! Now, the denominator for this ij
       do j=0,NY+1
         do k=0,NZM
           do i=0,NXM
             denominator(i,k,j)=Mij(i,k,j)*Mij(i,k,j)
           end do
         end do
       end do

        DO j=0,NY+1
          denominator_sum(j) = denominator_sum(j) +
     &                cross*SUM(denominator(0:NXM,0:NZM,j))
          numerator_sum(j)   = numerator_sum(j) +
     &                cross*SUM(numerator(0:NXM,0:NZM,j))
      ENDDO 

  
! End to ij
      end do



! Get plane average of numerator and denominator, add to C
! Note, since both the numerator and denominator are averaged, there
! is not need to divide the sum by NX*NZ

       do j=0,NY+1

          if (denominator_sum(j).ne.0.) then
            C_DYN(j) =-0.5d0* numerator_sum(j)/denominator_sum(j)


          else 
            C_DYN(j)=0.d0  
          end if
       end do




! We are now done with the dynamic procedure to calculate C
      

! If C_DYN < 0 at any level, set C_DYN=0 for numerical stability
      do j=0,NY+1
        if (C_DYN(j).lt.0) C_DYN(j)=0.d0
      end do

! At this point we have C_DYN and Sij
      

 
      
! End if compute C_DYN
       END IF

      

! Get the eddy viscosity at GY points
! NU_T = C_DYN * (DELTA^2)*|S|
! Use exact second order interpolation to get C_DYN and S1 interpolated to
! GY points
!      DO J=J1,J2
       DO J=J1i,min(NY+1,J2+1)
        DO K=0,NZM
          DO I=0,NXM
            NU_T(I,K,J)=
     &        (DYF(j-1)*C_DYN(j)+DYF(j)*C_DYN(j-1))/(2.d0*DY(j))
     &        *DELTA_Y(j)**2.d0
     &        *(DYF(j-1)*S1(i,k,j)+DYF(j)*S1(i,k,j-1))/(2.d0*DY(j))
          END DO
        END DO
      END DO
      
c      write(400,*) time_step, sum(NU_T(0:NXM,0:NZM,1)), 10000, J2
c      write(410,*) time_step, sum(NU_T(0:NXM,0:NZM,NY)), 20000, J2  


! Calculate TAUij in physical space, stored in Sij

      do ij=1,6
      if (LES_MODEL_TYPE.eq.2) then
! Dynamic Smagorinsky model, no scale similar part
      if ((ij.eq.1).or.(ij.eq.3).or.(ij.eq.5)) then
! Here, Sij is defined at GYF points 
      do j=0,NY+1
        do k=0,NZM
          do i=0,NXM
            Sij(i,k,j,ij)=
     &        -2.d0*C_DYN(j)*DELTA_YF(j)**2.d0*S1(i,k,j)*Sij(i,k,j,ij)
          end do 
        end do
      end do
c      do k=0,NZM
c        do i=0,NXM
! Get TAUij at ghost cells for computing derivative
c          Sij(i,k,0,ij)=Sij(i,k,1,ij)
c          Sij(i,k,NY+1,ij)=Sij(i,k,NY,ij)
c        end do
c      end do
      else if (ij.eq.2) then
! Sij(:,:,:,2) = du2/dy, but this term will be added implicitly through
! an eddy viscosity, so set it equal to zero here
        do j=0,NY+1
          do k=0,NZM
            do i=0,NXM
              Sij(i,k,j,ij)=0.d0
            end do
          end do
        end do 
      else if (ij.eq.4) then 
! Here, Sij is defined at GY points, interpolate C_DYN, etc 
! Use exact second order interpolation
! Sij(:,:,:,4)=0.5*(dU1/dy + dU2/dx)
! But dU1/dy will be accounted for as an implicit eddy viscosity term,
! So, subtract if off from Sij here
      do j=1,NY+1
        do k=0,NZM
          do i=0,NXM
            temp_1(i,k,j)=-2.d0
     &        *(DYF(j-1)*C_DYN(j)+DYF(j)*C_DYN(j-1))/(2.d0*DY(j))
     &        *DELTA_Y(j)**2.d0
     &        *(DYF(j-1)*S1(i,k,j)+DYF(j)*S1(i,k,j-1))/(2.d0*DY(j))
     &        *Sij(i,k,j,ij)
            Sij(i,k,j,ij)=-2.d0
     &        *(DYF(j-1)*C_DYN(j)+DYF(j)*C_DYN(j-1))/(2.d0*DY(j))
     &        *DELTA_Y(j)**2.d0
     &        *(DYF(j-1)*S1(i,k,j)+DYF(j)*S1(i,k,j-1))/(2.d0*DY(j))
     &        *(Sij(i,k,j,ij)-0.5d0*(U1(i,k,j)-U1(i,k,j-1))/DY(j))
          end do
        end do
      end do 
! Get TAUij at ghost cells for computing derivative
!      do k=0,NZM
!        do i=0,NXM
!          temp_1(i,k,1)   =temp_1(i,k,2)
!          temp_1(i,k,NY+1)=temp_1(i,k,NY)
!          Sij(i,k,1,ij)   =Sij(i,k,2,ij)
!          Sij(i,k,NY+1,ij)=Sij(i,k,NY,ij)
!        end do
!       end do
      else if (ij.eq.6) then     
! Here, Sij is defined at GY points, interpolate C_DYN, etc 
! Use exact second order interpolation
! Sij(:,:,:,6)=0.5*(dU3/dy + dU2/dz)
! But dU3/dy will be accounted for as an implicit eddy viscosity term,
! So, subtract if off from Sij here
      do j=1,NY+1
        do k=0,NZM
          do i=0,NXM
            temp_2(i,k,j)=-2.d0
     &        *(DYF(j-1)*C_DYN(j)+DYF(j)*C_DYN(j-1))/(2.d0*DY(j))
     &        *DELTA_Y(j)**2.d0
     &        *(DYF(j-1)*S1(i,k,j)+DYF(j)*S1(i,k,j-1))/(2.d0*DY(j))
     &        *Sij(i,k,j,ij)
            Sij(i,k,j,ij)=-2.d0
     &        *(DYF(j-1)*C_DYN(j)+DYF(j)*C_DYN(j-1))/(2.d0*DY(j))
     &        *DELTA_Y(j)**2.d0
     &        *(DYF(j-1)*S1(i,k,j)+DYF(j)*S1(i,k,j-1))/(2.d0*DY(j))
     &        *(Sij(i,k,j,ij)-0.5d0*(U3(i,k,j)-U3(i,k,j-1))/DY(j))
          end do
        end do
      end do 
! Get TAUij at ghost cells for computing derivative
!      do k=0,NZM
!        do i=0,NXM
!          temp_2(i,k,1)   =temp_2(i,k,2)
!          temp_2(i,k,NY+1)=temp_2(i,k,NY)
!          Sij(i,k,1,ij)=Sij(i,k,2,ij)
!          Sij(i,k,NY+1,ij)=Sij(i,k,NY,ij)
!        end do
!       end do

! End if ij
      end if

      else if (LES_MODEL_TYPE.eq.3) then
! Model type = 3, dynamic mixed model with scale similar part
! Always define temp at GYF points to match U_2BAR
! temp=Ui*Uj:
       SELECT CASE (ij)
       CASE(1)
         do j=0,NY+1
           do k=0,NZM
             do i=0,NXM 
               temp(i,k,j)=U1(i,k,j)*U1(i,k,j)
             end do
           end do
         end do
       CASE(2)
         do j=0,NY+1
           do k=0,NZM
             do i=0,NXM 
! Interpolate U2^2 to GYF points
               temp(i,k,j)=(0.5d0*(U2(i,k,j)+U2(i,k,j+1)))**2.d0
             end do
           end do
         end do
! Get interpolate temp to ghost cells 
!         do k=0,NZM
!           do i=0,NXM 
!             temp(i,k,0)=(U2(i,k,1)*(1.d0-DY(1)/(2.d0*DYF(1)))
!     &                   -U2(i,k,2)*DY(1)/(2.d0*DYF(1)))**2.d0
!             temp(i,k,NY+1)=(U2(i,k,NY+1)*(1.d0+DY(NY+1)
!     &        /(2.d0*DYF(NY)))-U2(i,k,NY)*DY(NY+1)/(2.d0*DYF(NY)))**2.d0
!           end do
!         end do
       CASE(3)
         do j=0,NY+1
           do k=0,NZM
             do i=0,NXM 
               temp(i,k,j)=U3(i,k,j)*U3(i,k,j)
             end do
           end do
         end do
       CASE(4) 
         do j=0,NY+1
           do k=0,NZM
             do i=0,NXM 
! Interpolate U2 to GYF points
               temp(i,k,j)=0.5d0*(U2(i,k,j)+U2(i,k,j+1))*U1(i,k,j)
             end do
           end do
         end do
! Get interpolate temp to ghost cells 
!         do k=0,NZM
!           do i=0,NXM 
!             temp(i,k,0)=(U2(i,k,1)*(1.d0-DY(1)/(2.d0*DYF(1)))
!     &                   -U2(i,k,2)*DY(1)/(2.d0*DYF(1)))*U1(i,k,0)
!             temp(i,k,NY+1)=(U2(i,k,NY+1)*(1.d0+DY(NY+1)
!     &        /(2.d0*DYF(NY)))-U2(i,k,NY)*DY(NY+1)/(2.d0*DYF(NY)))
!     &         *U1(i,k,NY+1)
!           end do
!         end do
       CASE(5)
         do j=0,NY+1
           do k=0,NZM
             do i=0,NXM 
               temp(i,k,j)=U1(i,k,j)*U3(i,k,j)
             end do
           end do
         end do
       CASE(6) 
         do j=0,NY+1
           do k=0,NZM
             do i=0,NXM 
! Interpolate U3 to GY points
               temp(i,k,j)=0.5d0*(U2(i,k,j)+U2(i,k,j+1))*U3(i,k,j)
             end do
           end do
         end do
! Get interpolate temp to ghost cells 
!         do k=0,NZM
!           do i=0,NXM 
!             temp(i,k,0)=(U2(i,k,1)*(1.d0-DY(1)/(2.d0*DYF(1)))
!     &                  -U2(i,k,2)*DY(1)/(2.d0*DYF(1)))*U3(i,k,0)
!             temp(i,k,NY+1)=(U2(i,k,NY+1)*(1.d0+DY(NY+1)
!     &        /(2.d0*DYF(NY)))-U2(i,k,NY)*DY(NY+1)/(2.d0*DYF(NY)))
!     &         *U3(i,k,NY+1)
!           end do
!         end do
       END SELECT

! Filter temp grid filter operation filter type = 1
       call les_filter_chan(temp,0,NY+1,1)

      if ((ij.eq.1).or.(ij.eq.3).or.(ij.eq.5)) then
! Here, Sij is defined at GYF points 

      do j=0,NY+1
        do k=0,NZM
          do i=0,NXM
            Sij(i,k,j,ij)=temp(i,k,j)
     &         -U_2BAR(i,k,j,U_index1(ij))*U_2BAR(i,k,j,U_index2(ij))
     &         -2.d0*C_DYN(j)*DELTA_YF(j)**2.d0*S1(i,k,j)*Sij(i,k,j,ij)
          end do
        end do 
      end do
! Get TAUij at ghost cells for computing derivative
!      do k=0,NZM
!        do i=0,NXM
!          Sij(i,k,0,ij)=Sij(i,k,1,ij)
!          Sij(i,k,NY+1,ij)=Sij(i,k,NY,ij)
!        end do
!       end do
       else if (ij.eq.2) then
! Sij(:,:,:,2) = du2/dy, but this term will be added implicitly through
! an eddy viscosity, so set the Smagorinsky pert equal to zero here
        do j=0,NY+1
          do k=0,NZM
            do i=0,NXM
              Sij(i,k,j,ij)=temp(i,k,j)
     &        -U_2BAR(i,k,j,U_index1(ij))*U_2BAR(i,k,j,U_index2(ij))
            end do
          end do
        end do 
      else if (ij.eq.4) then
! Here, Sij is defined at GY points, interpolate C_DYN, etc 
! Use exact second order interpolation
! Sij(:,:,:,4)=0.5*(dU1/dy + dU2/dx)
! But the dU1/dy term will be accounted for as an implicit eddy viscosity
! so subtract this term from Sij in the Smagorinsky part
      do j=1,NY+1
        do k=0,NZM
          do i=0,NXM
            temp_1(i,k,j)=
     &          (temp(i,k,j)*DYF(j-1)+temp(i,k,j-1)*DYF(j))/(2.d0*DY(j))
     &         -((U_2BAR(i,k,j,U_index1(ij))*DYF(j-1)
     &          +U_2BAR(i,k,j-1,U_index1(ij))*DYF(j))/(2.d0*DY(j)))
     &         *((U_2BAR(i,k,j,U_index2(ij))*DYF(j-1)
     &          +U_2BAR(i,k,j-1,U_index2(ij))*DYF(j))/(2.d0*DY(j)))
     &      -2.d0*(C_DYN(j)*DYF(j-1)+C_DYN(j-1)*DYF(j))/(2.d0*DY(j))
     &          *DELTA_Y(j)**2.d0
     &          *(S1(i,k,j)*DYF(j-1)+S1(i,k,j-1)*DYF(j))/(2.d0*DY(j))
     &          *Sij(i,k,j,ij)

            Sij(i,k,j,ij)=
     &         (temp(i,k,j)*DYF(j-1)+temp(i,k,j-1)*DYF(j))/(2.d0*DY(j))
     &         -((U_2BAR(i,k,j,U_index1(ij))*DYF(j-1)
     &          +U_2BAR(i,k,j-1,U_index1(ij))*DYF(j))/(2.d0*DY(j)))
     &         *((U_2BAR(i,k,j,U_index2(ij))*DYF(j-1)
     &          +U_2BAR(i,k,j-1,U_index2(ij))*DYF(j))/(2.d0*DY(j)))
     &      -2.d0*(C_DYN(j)*DYF(j-1)+C_DYN(j-1)*DYF(j))/(2.d0*DY(j))
     &          *DELTA_Y(j)**2.d0
     &          *(S1(i,k,j)*DYF(j-1)+S1(i,k,j-1)*DYF(j))/(2.d0*DY(j))
     &          *(Sij(i,k,j,ij)-0.5d0*(U1(i,k,j)-U1(i,k,j-1))/DY(j))
          end do
        end do
      end do 
! Get TAUij at ghost cells for computing derivative
!      do k=0,NZM
!        do i=0,NXM
!          temp_1(i,k,1)=temp_1(i,k,2)
!          temp_1(i,k,NY+1)=temp_1(i,k,NY)
!          Sij(i,k,1,ij)=Sij(i,k,2,ij)
!          Sij(i,k,NY+1,ij)=Sij(i,k,NY,ij)
!        end do
!       end do
      else if (ij.eq.6) then
! Here, Sij is defined at GY points, interpolate C_DYN, etc 
! Use exact second order interpolation
! Sij(:,:,:,6)=0.5*(dU3/dy + dU2/dz)
! But the dU1/dy term will be accounted for as an implicit eddy viscosity
! so subtract this term from Sij in the Smagorinsky part
      do j=1,NY+1
        do k=0,NZM
          do i=0,NXM
            temp_2(i,k,j)=
     &          (temp(i,k,j)*DYF(j-1)+temp(i,k,j-1)*DYF(j))/(2.d0*DY(j))
     &         -((U_2BAR(i,k,j,U_index1(ij))*DYF(j-1)
     &          +U_2BAR(i,k,j-1,U_index1(ij))*DYF(j))/(2.d0*DY(j)))
     &         *((U_2BAR(i,k,j,U_index2(ij))*DYF(j-1)
     &          +U_2BAR(i,k,j-1,U_index2(ij))*DYF(j))/(2.d0*DY(j)))
     &      -2.d0*(C_DYN(j)*DYF(j-1)+C_DYN(j-1)*DYF(j))/(2.d0*DY(j))
     &          *DELTA_Y(j)**2.d0
     &          *(S1(i,k,j)*DYF(j-1)+S1(i,k,j-1)*DYF(j))/(2.d0*DY(j))
     &          *Sij(i,k,j,ij)

            Sij(i,k,j,ij)=
     &         (temp(i,k,j)*DYF(j-1)+temp(i,k,j-1)*DYF(j))/(2.d0*DY(j))
     &         -((U_2BAR(i,k,j,U_index1(ij))*DYF(j-1)
     &          +U_2BAR(i,k,j-1,U_index1(ij))*DYF(j))/(2.d0*DY(j)))
     &         *((U_2BAR(i,k,j,U_index2(ij))*DYF(j-1)
     &          +U_2BAR(i,k,j-1,U_index2(ij))*DYF(j))/(2.d0*DY(j)))
     &      -2.d0*(C_DYN(j)*DYF(j-1)+C_DYN(j-1)*DYF(j))/(2.d0*DY(j))
     &          *DELTA_Y(j)**2.d0
     &          *(S1(i,k,j)*DYF(j-1)+S1(i,k,j-1)*DYF(j))/(2.d0*DY(j))
     &          *(Sij(i,k,j,ij)-0.5d0*(U3(i,k,j)-U3(i,k,j-1))/DY(j))
          end do
        end do
      end do 
! Get TAUij at ghost cells for computing derivative
!      do k=0,NZM
!        do i=0,NXM
!          temp_2(i,k,1)=temp_2(i,k,2)
!          temp_2(i,k,NY+1)=temp_2(i,k,NY)
!          Sij(i,k,1,ij)=Sij(i,k,2,ij)
!          Sij(i,k,NY+1,ij)=Sij(i,k,NY,ij)
!        end do
!       end do

! End if ij
       end if
! End if Mixed Model
       end if
! Convert TAUij, now stored in Sij to Fourier space
       call FFT_XZ_TO_FOURIER(Sij(0,0,0,ij),CSij(0,0,0,ij),0,NY+1)
! End do ij
       end do 
! Convert TAUij, now stored in Sij to Fourier space (for only u momentum equation)
       call FFT_XZ_TO_FOURIER(temp_1,ctemp_1,0,NY+1)
       call FFT_XZ_TO_FOURIER(temp_2,ctemp_2,0,NY+1) 
! End if LES_MODEL_TYPE dynamic Smagorinsky or Dynamic mixed model
      else
        stop
        write(6,*) 'Error, unsupported LES_MODEL_TYPE chosen'
      end if



      IF ((U_BC_YMIN.EQ.3).or.(W_BC_YMIN.EQ.3)) THEN
        IF ( J1i .EQ. 2) THEN    
           DO J =1,2
            DO K=0,TNKZ
             DO I=0,NKX
               DO ij=2,6,2
                CSij(I,K,J,ij)=0.d0
               ENDDO
                ctemp_1(I,K,J)=0.d0
                ctemp_2(I,K,J)=0.d0
             END DO
            END DO
           ENDDO
         ENDIF
      ENDIF

      IF ((U_BC_YMAX.EQ.3).or.(W_BC_YMAX.EQ.3)) THEN
        IF (J2e .EQ. NY )THEN
          DO J =NY,NY+1
           DO K=0,TNKZ
            DO I=0,NKX
             DO ij=2,6,2
              CSij(I,K,J,ij)=0.d0
             END DO
              ctemp_1(I,K,J)=0.d0
              ctemp_2(I,K,J)=0.d0
            END DO
           END DO
          ENDDO
        ENDIF
      ENDIF



! Zero Tau above J2+1 where it is not used
      
!      DO J=J2+2,NY+1
!        DO K=0,TNKZ
!          DO I=0,NKX
!            DO ij=1,6
c              CSij(I,K,J,ij)=0.d0
!            END DO
c              ctemp_1(I,K,J) =0.d0
c              ctemp_2(I,K,J) =0.d0 
!          END DO
!        END DO
!      END DO

! Now, add the subgrid scale forcing to CFi
! (This includes the subgrid scale stress as an explicit R-K term

c      write(400,*) time_step, CSij(0,0,1,1), 10000, J2
c      write(410,*) time_step, CSij(0,0,NY,1), 20000, J2     

      DO J=J1,J2
        DO K=0,TNKZ
          DO I=0,NKX
            CF1(I,K,J)=CF1(I,K,J)
     &                -CIKX(I)*CSij(I,K,J,1)
     &                -(CSij(I,K,J+1,4)-CSij(I,K,J,4))/DYF(j)
     &                -CIKZ(K)*CSij(I,K,J,5)
           CF3(I,K,J)=CF3(I,K,J)
     &                -CIKX(I)*CSij(I,K,J,5)
     &                -(CSij(I,K,J+1,6)-CSij(I,K,J,6))/DYF(J)
     &                -CIKZ(K)*CSij(I,K,J,3)   
          END DO
        END DO
      END DO
      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NKX
           CF2(I,K,J)=CF2(I,K,J)
     &                -CIKX(I)*ctemp_1(I,K,J)
     &                -(CSij(I,K,J,2)-CSij(I,K,J-1,2))/DY(j)
     &                -CIKZ(K)*ctemp_2(I,K,J)
          END DO
        END DO
      END DO

c         write(*,*) 'LES is running'
! Periodically, output mean quantities
      IF ((MOD(TIME_STEP,SAVE_STATS_INT).EQ.0).AND.(RK_STEP.EQ.1)) THEN
! Get plane averages

c        write(*,*) 'LES data are saving'
        do j=1,NY
          S1_mean(j)=0.d0
          NU_T_mean(j)=0.d0
          do i=0,NXM
          do k=0,NZM
            S1_mean(j)=S1_mean(j)+S1(I,K,J)
            NU_T_mean(j)=NU_T_mean(j)+NU_T(I,K,J)
          end do
          end do
          S1_mean(j)=S1_mean(j)/dble(NX*NZ)
          NU_T_mean(j)=NU_T_mean(j)/dble(NX*NZ)
        end do

c        open(42,file='mean_les.txt',form='formatted',status='unknown')

        open(42,file=FNAME,form='formatted',status='unknown',
     &  position='append')
        write(42,*) TIME_STEP,TIME,DELTA_T
        do j=1,NY
          write(42,420) j,GYF(J)
     &    ,(dble(CSij(0,0,J,ij)),ij=1,6),
     &       C_DYN(J)*DELTA_YF(J)**2.d0*S1_mean(J),
     &      NU_T_mean(J),C_DYN(J),
     &      C_DYN(J)*DELTA_YF(J)**2.d0 
        end do
       close(42)

420     format(I3,' ',11(F20.9,' '))

       do j=0,NY+1
        NU_T_MEAN(j)=SUM(NU_T(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do

!     Output SGS contribution to the TKE equation
!     LES term in the tke equation: -<u_i'dtau_ij'/dx_j>

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     sgs contribution to the transport: -dtau_ij'u_i'/dx_j 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!     i=1,j=1
!     First get Sij in physical space through temp

c      do j=2,NY
c        do k=0,TNKZ
c          do i=0,NKX
c            ctemp(i,k,j)=(DYF(J-1)*CSij(i,k,j,1)+DYF(J)*CSij(i,k,j-1,1))
c     &          /(2.d0*DY(j))
c          end do
c        end do
! Zero the mean mode to get tau'_ij 
c        ctemp(0,0,j)=0.d0
c      end do
c! Convert to Physical space
c      call fft_xz_to_physical(ctemp,temp,0,NY+1)
c! Multiply by U1'
c      do j=2,NY
c        do k=0,NZM
c          do i=0,NXM
c            temp(i,k,j)=temp(i,k,j)*(DYF(J-1)*(U1(i,k,j)-U1_bar(j))
c     &         +DYF(j)*(U1(i,k,j-1)-U1_bar(j-1)))/(2.d0*DY(J))
c          end do
c        end do
c      end do
c! Convert to Fourier space
c      call fft_xz_to_fourier(temp,ctemp,0,NY+1)
! Get x-derivative      
c      do j=2,NY
c        do k=0,TNKZ
c          do i=0,NKX
c            ctemp(i,k,j)=-CIKX(I)*ctemp(i,k,j)
c          end do
c        end do
c      end do
c! Convert to physical space
c      call FFT_XZ_TO_PHYSICAL(ctemp,temp,0,NY+1)
c
c      do j=2,NY
c        tke_sgs_t(j)=SUM(temp(0:NXM,0:NZM,j))/dble(NX*NZ)
c      end do      

! i=1,j=2
! First get Sij in physical space through temp
! define temp at GYF points
      do j=1,NY
        do k=0,TNKZ
          do i=0,NKX
c            ctemp(i,k,j)=0.5d0*(CSij(i,k,j,4)+CSij(i,k,j+1,4))
             ctemp(i,k,j)=0.5d0*(ctemp_1(i,k,j)+ctemp_1(i,k,j+1))
          end do
        end do
! Zero the mean mode to get tau'_ij 
        ctemp(0,0,j)=0.d0
      end do
! Convert to Physical space
      call fft_xz_to_physical(ctemp,temp,0,NY+1)
! Multiply by U1'
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            temp(i,k,j)=temp(i,k,j)*(U1(i,k,j)-U1_bar(j))
          end do
        end do
      end do
! Now, get a y-derivatives
      do j=NY,2,-1
        do k=0,NZM
          do i=0,NXM
            temp(i,k,j)=-(temp(i,k,j)-temp(i,k,j-1))/DY(j)
          end do
        end do
      end do
      do j=2,NY
        tke_sgs_t(j) =  SUM(temp(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do
c!i=1,j=3
c! First get Sij in physical space through temp
c      do j=2,NY
c        do k=0,TNKZ
c          do i=0,NKX
c            ctemp(i,k,j)= (DYF(J-1)*CSij(i,k,j,5)
c     &                  + DYF(J)*CSij(i,k,j-1,5))/(2.d0*DY(j))
c          end do
c        end do
c! Zero the mean mode to get tau'_ij
c        ctemp(0,0,j)=0.d0
c      end do
! Convert to Physical space
c      call fft_xz_to_physical(ctemp,temp,0,NY+1)
c! Multiply by U1'
c      do j=2,NY
c        do k=0,NZM
c          do i=0,NXM
c            temp(i,k,j)=temp(i,k,j)*(DYF(J-1)*(U1(i,k,j)-U1_bar(j))
c     &         +DYF(j)*(U1(i,k,j-1)-U1_bar(j-1)))/(2.d0*DY(J))
c          end do
c        end do
c      end do
c! Convert to Fourier space
c      call fft_xz_to_fourier(temp,ctemp,0,NY+1)
c! Get z-derivative
c      do j=2,NY
c        do k=0,TNKZ
c          do i=0,NKX
c            ctemp(i,k,j)= - CIKZ(K)*ctemp(i,k,j)
c          end do
c        end do
c      end do
c      call FFT_XZ_TO_PHYSICAL(ctemp,temp,0,NY+1)
c      do j=2,NY
c        tke_sgs_t(j)=tke_sgs_t(j)
c     &       +SUM(temp(0:NXM,0:NZM,j))/dble(NX*NZ)
c      end do
c! i=2, j=1
c! First get Sij in physical space through temp
c      do j=2,NY
c        do k=0,TNKZ
c          do i=0,NKX
c            ctemp(i,k,j)=CSij(i,k,j,4)
c          end do
c        end do
c! Zero the mean mode to get tau'_ij
c        ctemp(0,0,j)=0.d0    
c      end do
c! Convert to Physical space
c      call fft_xz_to_physical(ctemp,temp,0,NY+1)
c! Multiply by U2'
c      do j=2,NY
c        do k=0,NZM
c          do i=0,NXM
c            temp(i,k,j)=temp(i,k,j)*(U2(i,k,j)-U2_bar(j))
c          end do
c        end do
c      end do
c! Convert to Fourier space
c      call fft_xz_to_fourier(temp,ctemp,0,NY+1)
c! Get x-derivative
c      do j=2,NY
c        do k=0,TNKZ
c          do i=0,NKX
c            ctemp(i,k,j)= - CIKX(I)*ctemp(i,k,j)
c          end do
c        end do
c      end do
c      call FFT_XZ_TO_PHYSICAL(ctemp,temp,0,NY+1)
c      do j=2,NY
c        tke_sgs_t(j)=tke_sgs_t(j)
c    &         +SUM(temp(0:NXM,0:NZM,j))/dble(NX*NZ)
c      end do

! i=2,j=2
! First get Sij in physical space through temp
! First, define temp at GYF points
      do j=1,NY
        do k=0,TNKZ
          do i=0,NKX
            ctemp(i,k,j)=CSij(i,k,j,2)
          end do
       end do
      enddo
      call fft_xz_to_physical(ctemp,temp,0,NY+1)
        do j=1,NY
         do k=0,NZM
          do i=0,NXM
            temp(i,k,j) = temp(i,k,j)-2.d0*C_DYN(j)*DELTA_YF(j)**2.d0
     &                      *S1(i,k,j)*st_rij(i,k,j,2)
          enddo
         enddo
        enddo
       call fft_xz_to_fourier(temp,ctemp,0,NY+1)
      do j=1,NY
! Zero the mean mode to get tau'_ij
        ctemp(0,0,j)=0.d0
      end do
! Convert to Physical space
      call fft_xz_to_physical(ctemp,temp,0,NY+1)
! Multiply by U2'
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            temp(i,k,j)=temp(i,k,j)*0.5*(U2(i,k,j)+U2(i,k,j+1)
     &                   -U2_bar(J)-U2_bar(J+1))
          end do
        end do
      end do
! Take vertical derivative of temp, re-defining at GY points
      do j=NY,2,-1
        do k=0,NZM
          do i=0,NXM
            temp(i,k,j)=-(temp(i,k,j)-temp(i,k,j-1))/DY(j)
          end do
        end do
      end do
      do j=2,NY
        tke_sgs_t(j)=tke_sgs_t(j)
     &           +SUM(temp(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do
c! i=2,j=3
c! First get Sij in physical space through temp
c      do j=2,NY
c        do k=0,TNKZ
c          do i=0,NKX
c            ctemp(i,k,j)=CSij(i,k,j,6)
c          end do
c        end do
c! Zero the mean mode to get tau'_ij
c        ctemp(0,0,j)=0.d0
c      end do
! Convert to Physical space
c      call fft_xz_to_physical(ctemp,temp,0,NY+1)
c! Multiply by U2'
c      do j=2,NY
c        do k=0,NZM
c          do i=0,NXM
c            temp(i,k,j)=temp(i,k,j)*(U2(i,k,j)-U2_bar(j))
c          end do
c        end do
c      end do
c! Convert to Fourier space
c      call fft_xz_to_fourier(temp,ctemp,0,NY+1)
c! Get z-derivative
c      do j=2,NY
c        do k=0,TNKZ
c          do i=0,NKX
c            ctemp(i,k,j)=-CIKZ(K)*ctemp(i,k,j)
c          end do
c        end do
c      end do
c      call FFT_XZ_TO_PHYSICAL(ctemp,temp,0,NY+1)
c      do j=2,NY
c        tke_sgs_t(j)=tke_sgs_t(j)
c     &       +SUM(temp(0:NXM,0:NZM,j))/dble(NX*NZ)
c      end do
c! i=3, j=1
c! First get Sij in physical space through temp
c      do j=2,NY
c        do k=0,TNKZ
c          do i=0,NKX
c            ctemp(i,k,j)=(DYF(j-1)*CSij(i,k,j,5)
c     &                   +DYF(J)*CSij(i,k,j-1,5))/(2.d0*DY(j))
c          end do
c        end do
c! Zero the mean mode to get tau'_ij
c       ctemp(0,0,j)=0.d0
c      end do
c! Convert to Physical space
c      call fft_xz_to_physical(ctemp,temp,0,NY+1)
c! Multiply by U3'
c      do j=2,NY
c        do k=0,NZM
c          do i=0,NXM
c            temp(i,k,j)=temp(i,k,j)*(DYF(J-1)*(U3(i,k,j)-U3_bar(j))
c     &         +DYF(j)*(U3(i,k,j-1)-U3_bar(j-1)))/(2.d0*DY(J))
c          end do
c        end do
c      end do
c! Convert to Fourier space
c      call fft_xz_to_fourier(temp,ctemp,0,NY+1)
c! Get x-derivative
c      do j=2,NY
c        do k=0,TNKZ
c          do i=0,NKX
c            ctemp(i,k,j)= - CIKX(I)*ctemp(i,k,j)
c          end do
c        end do
c      end do
c      call FFT_XZ_TO_PHYSICAL(ctemp,temp,0,NY+1)
c      do j=2,NY
c        tke_sgs_t(j) = tke_sgs_t(j)
c     &        +SUM(temp(0:NXM,0:NZM,j))/dble(NX*NZ)
c      end do
! i=3, j=2
! First get Sij in physical space through temp
! First get temp at GYF points
      do j=1,NY
        do k=0,TNKZ
          do i=0,NKX
            ctemp(i,k,j)=0.5d0*(ctemp_2(i,k,j)+ctemp_2(i,k,j+1))
          end do
        end do
! Zero the mean mode to get tau'_ij
        ctemp(0,0,j)=0.d0
      end do
! Convert to Physical space
      call fft_xz_to_physical(ctemp,temp,0,NY+1)
! Multiply by U3'
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            temp(i,k,j)=temp(i,k,j)*(U3(i,k,j)-U3_bar(j))
          end do
        end do
      end do
! Now, re-define temp at GY points
      do j=NY,2,-1
        do k=0,NZM
          do i=0,NXM
            temp(i,k,j)=-(temp(i,k,j)-temp(i,k,j-1))/DY(j)
          end do
        end do
      end do
      do j=2,NY
        tke_sgs_t(j)=tke_sgs_t(j)
     &        + SUM(temp(0:NXM,0:NZM,j))/dble(NX*NZ)
      end do
c! i=3, j=3
c! First get Sij in physical space through temp
c      do j=2,NY
c        do k=0,TNKZ
c          do i=0,NKX
c            ctemp(i,k,j)=(DYF(j-1)*CSij(i,k,j,3)
c     &                   +DYF(J)*CSij(i,k,j-1,3))/(2.d0*DY(j))
c          end do
c        end do
c! Zero the mean mode to get tau'_ij
c        ctemp(0,0,j)=0.d0           
c      end do
c! Convert to Physical space
c      call fft_xz_to_physical(ctemp,temp,0,NY+1)
c! Multiply by U3'
c      do j=2,NY
c        do k=0,NZM
c          do i=0,NXM
c            temp(i,k,j)=temp(i,k,j)*(DYF(J-1)*(U3(i,k,j)-U3_bar(j))
c     &         +DYF(j)*(U3(i,k,j-1)-U3_bar(j-1)))/(2.d0*DY(J))
c          end do
c        end do
c      end do
c! Convert to Fourier space
c      call fft_xz_to_fourier(temp,ctemp,0,NY+1)
c! Get z-derivative
c      do j=2,NY
c        do k=0,TNKZ
c          do i=0,NKX
c            ctemp(i,k,j)=-CIKZ(K)*ctemp(i,k,j)
c          end do
c        end do
c      end do
c      call FFT_XZ_TO_PHYSICAL(ctemp,temp,0,NY+1)
c      do j=2,NY
c        tke_sgs_t(j)=tke_sgs_t(j)
c     &        + SUM(temp(0:NXM,0:NZM,j))/dble(NX*NZ)
c      end do
c

     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     sgs contribution to the production: -<tau_ij><S_ij>
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      

!     To get S_ij for strain rate tensor we need to store tau_ij 
  
      DO ij = 1,6
      
        do j=1,NY
         do k=0,TNKZ
          do i=0,NKX
            ctemp(i,k,j) = CSij(i,k,j,ij)   
          enddo
         enddo
        enddo

      if ( (ij.eq. 1).or.(ij.eq.3).or.(ij.eq.5)) then
        call fft_xz_to_physical(ctemp,temp,0,NY+1)
        do j=1,NY
         do k=0,NZM
          do i=0,NXM
            Sij(i,k,j,ij) = temp(i,k,j)
          enddo
         enddo
        enddo
      else if (ij.eq.2) then
! we need to take into account Smagorinsky part du2/dy
        call fft_xz_to_physical(ctemp,temp,0,NY+1)
        do j=1,NY
         do k=0,NZM
          do i=0,NXM
            Sij(i,k,j,ij) = temp(i,k,j)-2.d0*C_DYN(j)*DELTA_YF(j)**2.d0
     &                      *S1(i,k,j)*st_rij(i,k,j,ij)
          enddo
         enddo
        enddo
      else if (ij.eq.4) then
        call fft_xz_to_physical(ctemp_1,temp_1,0,NY+1)
        do j=1,NY
         do k=0,NZM
          do i=0,NXM
            Sij(i,k,j,ij) = temp_1(i,k,j)
          enddo
         enddo
        enddo
      else if (ij.eq.6) then
        call fft_xz_to_physical(ctemp_2,temp_2,0,NY+1)
        do j=1,NY
         do k=0,NZM
          do i=0,NXM
            Sij(i,k,j,ij) = temp_2(i,k,j)
          enddo
         enddo
        enddo
      endif
       
      ENDDO
      
! Need to calculate Strain Rate

      do j=0,NY+1
        do ij=1,6
         Sij_mean(j,ij) = SUM(St_rij(0:NXM,0:NZM,j,ij))/dble(NX*NZ)
         TAU_mean(j,ij) = SUM(Sij(0:NXM,0:NZM,j,ij))/dble(NX*NZ)
        enddo
      enddo
         
c      do j=1,NY
c          check_tau(j) = 0.d0
c          check_tau(j)=SUM(temp_1(0:NXM,0:NZM,j))/dble(NX*NZ)
c      enddo 



      do j=2,NY 
         tke_sgs_p(j) = 0.d0
         do ij=1,3 
       tke_sgs_p(j)=tke_sgs_p(j)-(DYF(J-1)*Sij_mean(j,ij)*TAU_mean(j,ij)
     &         +DYF(j)*Sij_mean(j-1,ij)*TAU_mean(j-1,ij))/(2.d0*DY(J))
         enddo 
      enddo    

      
      do j=2,NY
       tke_sgs_p(j)= tke_sgs_p(j)- 2.0*Sij_mean(j,4)*TAU_mean(j,4)  
     &                           - 2.0*Sij_mean(j,6)*TAU_mean(j,6)
     &                 - 2.0*(DYF(J-1)*Sij_mean(j,5)*TAU_mean(j,5)
     &         +DYF(j)*Sij_mean(j-1,5)*TAU_mean(j-1,5))/(2.d0*DY(J))  
      enddo
        
         


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     sgs contribution to the dissipation: -<tau_ij*S_ij>
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      do j=2,NY
         tke_sgs_diss(j) = 0.d0
         do ij=1,6
          if ( (ij .eq. 4) .or. (ij .eq. 6) ) then
           do k=0,NZM
            do i=0,NXM
             temp(i,k,j) =   2.0*Sij(i,k,j,ij)*st_rij(i,k,j,ij)
            enddo
           enddo
          elseif ( ij .eq. 5 ) then
           do k=0,NZM
            do i=0,NXM
             temp(i,k,j)=  2.0*(DYF(J-1)*Sij(i,k,j,ij)*st_rij(i,k,j,ij)
     &     + DYF(j)*Sij(i,k,j-1,ij)*st_rij(i,k,j-1,ij))/(2.d0*DY(J))
            enddo
           enddo
          else
            do k=0,NZM
            do i=0,NXM
             temp(i,k,j)=  (DYF(J-1)*Sij(i,k,j,ij)*st_rij(i,k,j,ij)
     &     + DYF(j)*Sij(i,k,j-1,ij)*st_rij(i,k,j-1,ij))/(2.d0*DY(J))
            enddo
           enddo
          endif
          
         
       tke_sgs_diss(j)=tke_sgs_diss(j)+SUM(temp(0:NXM,0:NZM,j))
     &            /dble(NZ*NX) 
       enddo

       enddo      


 
     

!     Write tke_sgs component to the file
       open(52,file=FNAME_TKE,form='formatted',status='unknown',
     &  position='append')

       write(52,*)TIME_STEP,TIME,DELTA_T
        do j=1,NY
         write(52,520)J,GY(J),tke_sgs_t(J),tke_sgs_p(J),
     &                tke_sgs_diss(J)
        enddo
       
       close(52)
     
520   format(I3,'',4(F20.9,''))
      END IF


      RETURN
      END



      subroutine compute_strain_chan
C This subroutine computes S_ij for the filtered velocity field
C The input velocity field should be in fourier space in the periodic
C directions.
C For use in the LES model in channel flow (2 periodic directions)
      include 'header_les'

      integer I,J,K,ij

       
      DO J=0,NY+1
        DO K=0,TNKZ
          DO I=0,NKX
            CSij(I,K,J,1)=CIKX(I)*CU1(I,K,J)
            CSij(I,K,J,2)=(CU2(I,K,J+1)-CU2(I,K,J))/DYF(j) 
            CSij(I,K,J,3)=CIKZ(K)*CU3(I,K,J)  
            CSij(I,K,J,5)=0.5d0*(CIKZ(K)*CU1(I,K,J)+CIKX(I)*CU3(I,K,J))
          END DO
        END DO
      END DO
      DO J=0,NY+2
        DO K=0,TNKZ
          DO I=0,NKX
            CSij(I,K,J,4)=0.5d0*( (CU1(I,K,J)-CU1(I,K,J-1))/DY(j)
     &                          +CIKX(I)*CU2(I,K,J) ) 
            CSij(I,K,J,6)=0.5d0*( CIKZ(K)*CU2(I,K,J)
     &                         +(CU3(I,K,J)-CU3(I,K,J-1))/DY(j) )
          END DO
        END DO
      END DO  


! Convert rate of strain tensor to physical space
      do ij=1,6
        call FFT_XZ_TO_PHYSICAL(CSij(0,0,0,ij),Sij(0,0,0,ij),0,NY)
        call FFT_XZ_TO_PHYSICAL(CSij(0,0,NY+1,ij),Sij(0,0,NY+1,ij),0,1)
      end do 
! We now have S_ij in Physical space
      RETURN
      END

 
      subroutine les_filter_chan(A,jstart,jend, filter_type)
! This subroutine applies the les filter to the input field
! The indices to the start and end of the array in the y-direction
! are also inputted to make the routine cablable of filtering fields
! at either GYF or GY points.
! The array that is passed should be in physical space
      integer jstart,jend,NX,NY,NZ,NXM,NZM,N_TH,filter_type
      INCLUDE 'grid_def'

      real*8 A(0:NX+1,0:NZ+1,0:NY+1)
      real*8 B(0:NX-1,0:NZ-1,0:NY+1)

      integer im2(0:NX-1),im1(0:NX-1),ip1(0:NX+1),ip2(0:NX+2)
      integer km2(0:NZ-1),km1(0:NZ-1),kp1(0:NZ+1),kp2(0:NZ+2)

! These are the weights for the filtering operation used
      real*8 W0,W1,W2,Wm1,Wm2,Wm1_j,W0_j,W1_j

! filter type = 1: grid filter operation
! filter type = 2: test filter operation
! The following is for the 3-point trapezoidal rule, alpha*beta=sqrt(6)
      
      If ( filter_type .eq. 1 ) then
       Wm2=0.d0
       Wm1=1.d0/8.d0
       W0=6.d0/8.d0
       W1=1.d0/8.d0
       W2=0.d0
      elseif ( filter_type .eq. 2 ) then
!     else 
!      Wm1_j=1.d0/4.d0  
!      W0_j=1.d0/2.d0
!      W1_j=1.d0/4.d0
! The following is for the 5-point trapezoidal rule, alpha*beta=9
!       Wm2=1.d0/8.d0
!       Wm1=1.d0/4.d0
!       W0=1.d0/4.d0
!       W1=1.d0/4.d0
!       W2=1.d0/8.d0

       Wm2=0.d0
       Wm1=1.d0/4.d0
       W0=1.d0/2.d0
       W1=1.d0/4.d0
       W2=0.d0 
      else
       stop
       write(6,*) 'Error, unsupported LES_FILTER_TYPE chosen'
      endif

      NXM=NX-1
      NZM=NZ-1

!      do j=0,NY+1
!        do k=0,NZM
!          do i=0,NXM
!            B(i,k,j)=A(i,k,j)
!          end do
!        end do
!      end do

! Filter in the periodic directions using cshift
! Note, cshift if not used since it appears to be slower
! Apply filter in the x-direction
!      B=Wm2*CSHIFT(B,-2,1)+Wm1*CSHIFT(B,-1,1)+W0*B+W1*CSHIFT(B,1,1)
!     &       +W2*CSHIFT(B,2,1)

! Filter using more efficient F77 syntax:
! Set up array to loop around periodic directions
      do i=2,NXM
        im2(i)=i-2
      end do
      im2(1)=NXM
      im2(0)=NX-2
      do i=1,NXM
        im1(i)=i-1
      end do
      im1(0)=NXM
      do i=0,NX-2
        ip1(i)=i+1
      end do
      ip1(NXM)=0
      do i=0,NX-3
        ip2(i)=i+2    
      end do
      ip2(NX-2)=0
      ip2(NXM)=1

      do j=jstart,jend
        do k=0,NZM
          do i=0,NXM
            B(i,k,j)=Wm2*A(im2(i),k,j)+Wm1*A(im1(i),k,j)+W0*A(i,k,j)
     &         +W1*A(ip1(i),k,j)+W2*A(ip2(i),k,j)
          end do
        end do  
      end do
 
! Apply filter in the z-diretion
!      B=Wm2*CSHIFT(B,-2,2)+Wm1*CSHIFT(B,-1,2)+W0*B+W1*CSHIFT(B,1,2)
!     &       +W2*CSHIFT(B,2,2)
! Filter using more efficient F77 syntax:
! Set up array to loop around periodic directions
      do k=2,NZM
        km2(k)=k-2
      end do
      km2(1)=NZM
      km2(0)=NZ-2
      do k=1,NZM
        km1(k)=k-1
      end do
      km1(0)=NZM
      do k=0,NZ-2
        kp1(k)=k+1
      end do
      kp1(NZM)=0
      do k=0,NZ-3
        kp2(k)=k+2    
      end do
      kp2(NZ-2)=0
      kp2(NZM)=1

      do j=jstart,jend
        do k=0,NZM
          do i=0,NXM
            A(i,k,j)=Wm2*B(i,km2(k),j)+Wm1*B(i,km1(k),j)+W0*B(i,k,j)
     &         +W1*B(i,kp1(k),j)+W2*B(i,kp2(k),j)
          end do
        end do  
      end do

! Apply filter in the vertical direction at all physical cells
! (filter is not applied to ghost cells, but the values of the ghost cells
! is used for averaging)
!      B(:,:,jstart+1:jend-1)=W0*B(:,:,jstart:jend-2)
!     &                      +W1*B(:,:,jstart+1:jend-1)
!     &                      +W2*B(:,:,jstart+2:jend)
! Use more efficient F77 syntax:
!       do j=jstart+1,jend-1
!         do k=0,NZM
!           do i=0,NXM
!             B(i,k,j)=Wm1_j*B(i,k,j-1)+W0_j*B(i,k,j)+W1_j*B(i,k,j+1)  
!           end do
!         end do
!       end do

!      do j=jstart,jend
!        do k=0,NZM
!          do i=0,NXM
!            A(i,k,j)=B(i,k,j)
!          end do
!        end do
!      end do

      return
      end



      subroutine les_filter_chan_fourier(A,jstart,jend)
! This subroutine applies the les filter to the input field
! The filter is a spectral cutoff filter
! The indices to the start and end of the array in the y-direction
! are also inputted to make the routine cablable of filtering fields
! at either GYF or GY points.
! The array that is passed should be in physical space
      integer jstart,jend,NX,NZ,NY,NXM,NZM,i,j,k,N_TH

      INCLUDE 'grid_def'

      real*8 PI
      integer NKX,NKZ,TNKZ

      real*8 KX(0:NX/3),KZ(0:2*(NZ/3))  

      real*8 A(0:NX+1,0:NZ+1,0:NY+1)

      real alpha

      real*8 B(0:NX+1,0:NZ+1,0:NY+1)

      complex*16 CB(0:NX/2,0:NZ+1,0:NY+1)

      equivalence (B,CB)

! Set the ratio of filter scales
      parameter (alpha=2.d0)

      NXM=NX-1
      NZM=NZ-1


      PI = 4. * ATAN(1.0)

      LX=PI
      LZ=2.d0*PI
   
! Get the wavenumber vectors:
        NKX=NX/3
        DO I=0,NKX
          KX(I)=I*(2.*PI)/LX
        END DO

        NKZ=NZ/3
        TNKZ=NKZ*2
        DO K=0,NKZ
          KZ(K)=K*(2.*PI)/LZ
        END DO
        DO K=1,NKZ
          KZ(TNKZ+1-K)=-K*(2.*PI)/LZ
        END DO

       do j=0,NY+1
         do k=0,NZM
           do i=0,NXM
             B(i,k,j)=A(i,k,j)
           end do
         enddo
       end do 


! Convert to fourier space
      call fft_xz_to_fourier(B,CB,jstart,jend)

! Perform the filtering
      do j=jstart,jend
        do k=0,TNKZ
          do i=0,NKX
            if (sqrt(KX(I)**2.d0+KZ(K)**2.d0)
     &       .gt.sqrt(KX(NKX)**2.d0+KZ(NKZ)**2.d0)/alpha) then
                CB(i,k,j)=0.d0
            end if
          end do
        end do
      end do

! Now, convert back to physical space
      call fft_xz_to_physical(CB,B,jstart,jend)

       do j=jstart,jend 
         do k=0,NZM
           do i=0,NXM
             A(i,k,j)=B(i,k,j)
           end do
         end do
       end do 

      return
      end





