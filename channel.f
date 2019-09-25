C******************************************************************************|
C channel.f, the channel-flow solvers for zodiac.                  VERSION 0.9
C******************************************************************************|
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE INIT_CHAN
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Initialize any constants here
      INCLUDE 'header'
      INTEGER J, N
!
! Open a NetCDF file for writing standardized output

      PI=4.D0*ATAN(1.D0)

        write(*,*) 'U_BC_YMIN: ',U_BC_YMIN
        IF (U_BC_YMIN.EQ.0) THEN
          JSTART=2
        ELSE IF (U_BC_YMIN.EQ.1) THEN
          JSTART=1
        ELSE
          JSTART=2
        END IF
! Now, set the indexing for the scalar equations
        DO N=1,N_TH
          IF (TH_BC_YMIN(N).EQ.0) THEN
            JSTART_TH(N)=2
          ELSE IF (TH_BC_YMIN(N).EQ.1) THEN
            JSTART_TH(N)=0
          ELSE IF (TH_BC_YMIN(N).EQ.3) THEN
            JSTART_TH(N)=0
          ELSE
            JSTART_TH(N)=2
          END IF
        END DO
        write(*,*) 'U_BC_YMAX: ',U_BC_YMAX
        IF (U_BC_YMAX.EQ.0) THEN
          JEND=NY-1
        ELSE IF (U_BC_YMAX.EQ.1) THEN
          JEND=NY
        ELSE
          JEND=NY-1
        END IF

! Set the upper and lower limits of timestepping of the scalar equations
        DO N=1,N_TH
        IF (TH_BC_YMAX(N).EQ.0) THEN
          JEND_TH(N)=NY-1
        ELSE IF (TH_BC_YMAX(N).EQ.1) THEN
          JEND_TH(N)=NY+1
        ELSE
          JEND_TH(N)=NY-1
        END IF
        END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE RK_CHAN_1
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Main time-stepping algorithm for the channel-flow case.
C This algorithm uses Crank-Nicolson for all terms involving vertical
C derivatives (viscous and nonlinear) and 3rd order Runge-Kutta for the
C rest of the terms
C INPUTS  (in Fourier space):  CUi, CP, and (if k>1) CFi at (k-1)  (for i=1,2,3)
C OUTPUTS (in Fourier space):  CUi, CP, and (if k<3) CFi at (k)
C Each RK step, there are 14 FFT calls. 11 storage variables are used.     
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      INTEGER I,J,K,N      
      REAL*8 TEMP1, TEMP2, TEMP3, TEMP4, TEMP5, UBULK


      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE RK_CHAN_2
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Alternative time-stepping algorithm for the channel-flow case.
C This algorithm uses Crank-Nicolson for all viscous terms and 
C third order Runge-Kutta for all nonlinear terms
C INPUTS  (in Fourier space):  CUi, P, and (if k>1) CFi at (k-1)  (for i=1,2,3)
C OUTPUTS (in Fourier space):  CUi, P, and (if k<3) CFi at (k)
C Each RK step, there are 11 FFT calls. 11 storage variables are used.     
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      INCLUDE 'header'

      INTEGER I,J,K,N      
      REAL*8 TEMP1, TEMP2, TEMP3, TEMP4, TEMP5, UBULK

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE RK_CHAN_3
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Main time-stepping algorithm for the channel-flow case.
C This algorithm uses Crank-Nicolson for all vertical viscous terms
C derivatives (viscous and nonlinear) and 3rd order Runge-Kutta for the
C rest of the terms
C INPUTS  (in Fourier space):  CUi, CP, and (if k>1) CFi at (k-1)  (for i=1,2,3)
C OUTPUTS (in Fourier space):  CUi, CP, and (if k<3) CFi at (k)
C Each RK step, there are 14 FFT calls. 11 storage variables are used.     
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'
      INTEGER I,J,K,N      
      REAL*8 TEMP1, TEMP2, TEMP3, TEMP4, TEMP5, 
     &       TEMP2_v,  UBULK
      
      

C Communicate the information between ghost cells



 

      CALL GHOST_CHAN_MPI

C Define the constants that are used in the time-stepping
C For reference, see Numerical Renaissance
      TEMP1=NU * H_BAR(RK_STEP) / 2.0
      TEMP2=H_BAR(RK_STEP) / 2.0
      TEMP2_v=H_BAR(RK_STEP)
      TEMP3=ZETA_BAR(RK_STEP) * H_BAR(RK_STEP)
      TEMP4=H_BAR(RK_STEP)
      TEMP5=BETA_BAR(RK_STEP) * H_BAR(RK_STEP)

C First, we will compute the explicit RHS terms and store in Ri
C Note, Momentum equation and hence the RHS is evaluated at the
C corresponding velocity points.

C Store the old velocity in the RHS vector
      DO J=JSTART,JEND
        DO K=0,TNKZ
          DO I=0,NKX
            CR1(I,K,J)=CU1(I,K,J)
            CR3(I,K,J)=CU3(I,K,J)
          END DO
        END DO
      END DO
      DO J=2,NY 
        DO K=0,TNKZ
          DO I=0,NKX
            CR2(I,K,J)=CU2(I,K,J)
          END DO
        END DO
      END DO

C Add the R-K term from the rk-1 step 
      IF (RK_STEP .GT. 1) THEN
        DO J=JSTART,JEND
          DO K=0,TNKZ
            DO I=0,NKX
              CR1(I,K,J)=CR1(I,K,J)+TEMP3*CF1(I,K,J)
              CR3(I,K,J)=CR3(I,K,J)+TEMP3*CF3(I,K,J)
            END DO
          END DO
        END DO
        DO J=2,NY
          DO K=0,TNKZ
            DO I=0,NKX
              CR2(I,K,J)=CR2(I,K,J)+TEMP3*CF2(I,K,J)
            END DO
          END DO
        END DO
      END IF
          
      DO J=0,NY+1
        DO K=0,TNKZ
          DO I=0,NKX
            CF1(I,K,J)=0.d0
            CF2(I,K,J)=0.d0
            CF3(I,K,J)=0.d0
           END DO
         END DO
       END DO




C Take the y-derivative of the pressure at GY points in Fourier space
      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NKX
            CS1(I,K,J)=(CP(I,K,J) - CP(I,K,J-1)) / DY(J)
          END DO
        END DO
      END DO

C Add the pressure gradient to the RHS as explicit Euler
      DO J=JSTART,JEND
        DO K=0,TNKZ
          DO I=0,NKX
            CR1(I,K,J)=CR1(I,K,J)-TEMP4*(CIKX(I)*CP(I,K,J))
            CR3(I,K,J)=CR3(I,K,J)-TEMP4*(CIKZ(K)*CP(I,K,J))
          END DO
        END DO
      END DO
      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NKX
            CR2(I,K,J)=CR2(I,K,J)-TEMP4*CS1(I,K,J)
          END DO
        END DO
      END DO

C Here, add the constant, forcing pressure gradient
C There are several ways of doing this
      IF (F_TYPE.EQ.1) THEN 
C Add forcing for a constant pressure gradient
        DO J=JSTART,JEND
          CR1(0,0,J)=CR1(0,0,J)-TEMP4*PX0
        END DO
!        write(*,*) TEMP4, PX0
      ELSE IF (F_TYPE.EQ.0) THEN  
C Add the mean pressure gradient to keep Ubulk constant
C This section needs to be parallelized 
        DO J=JSTART,JEND
          CR3(0,0,J)=CR3(0,0,J)-TEMP4*PX0
        END DO

      ELSE IF (F_TYPE.EQ.2) THEN
C If oscillatory pressure gradient
        DO J=JSTART,JEND
          CR1(0,0,J)=CR1(0,0,J)-TEMP4*(PX0+AMP_OMEGA0*cos(OMEGA0*TIME))
        END DO
      
      ELSE IF (F_TYPE.EQ.4) THEN
C If oscillatory pressure gradient
        DO J=JSTART,JEND
          CR1(0,0,J)=CR1(0,0,J)-TEMP4*(PX0+AMP_OMEGA0*cos(OMEGA0*TIME))*
     :           cos(ANG_BETA)*exp(-(b_fact*(GYF(J)-GYF(JSTART)))**2.0)
        END DO 

C End if forcing type
      END IF
 
C Now compute the term R-K term Ai
C Compile terms of Ai in CFi which will be saved for next time step
C First, store the horizontal viscous terms in CFi
      DO J=JSTART,JEND
        DO K=0,TNKZ
          DO I=0,NKX
            CF1(I,K,J)=-NU * KX2(I) * CU1(I,K,J) 
     &            - NU * KZ2(K) * CU1(I,K,J)
            CF3(I,K,J)=-NU * KX2(I) * CU3(I,K,J) 
     &            - NU * KZ2(K) * CU3(I,K,J)
          END DO
        END DO
      END DO
      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NKX
            CF2(I,K,J)=-NU * KX2(I) * CU2(I,K,J) 
     &            - NU * KZ2(K) * CU2(I,K,J)
          END DO 
        END DO
      END DO

!   Add sponge in for open boundary condition

c       IF (J2e .eq. NY) THEN
c         CALL sponge    
c       ENDIF/Ratio_gr

!      write(6,*) 'I am here 2'
 
! Do for each scalar
      DO N=1,N_TH

      IF (F_TYPE .EQ. 4) THEN
        DO J=2,NY
         DO K=0,TNKZ
           DO I=0,NKX
! Use second order interpolation
              CF2(I,K,J)=CF2(I,K,J) - Ratio_gr*cos(ANG_BETA)*
     &       (CTH(I,K,J,N)*DYF(J-1)+CTH(I,K,J-1,N)*DYF(J))/(2.d0*DY(J))
           END DO
         END DO
        END DO

        DO J=JSTART,JEND
         DO K=0,TNKZ
          DO I=0,NKX
             CF1(I,K,J)=CF1(I,K,J)-Ratio_gr*sin(ANG_BETA)*CTH(I,K,J,N)
          END DO
         END DO
        END DO
      ELSE 
! If a scalar contributes to the denisty, RI_TAU is not equal to zero and
! add the buoyancy term as explicit R-K.  Don't add the 0,0 mode which 
! corresponds to a plane average.  The plane averaged density balances
! the hydrostratic pressure component.
      IF (DEV_BACK_TH) THEN

      IF(N_TH == 1) THEN
       DO J=2,NY
         DO K=0,TNKZ
           DO I=0,NKX
! Use second order interpolation
              CF2(I,K,J)=CF2(I,K,J) + Ratio_gr*
     &       (CTH(I,K,J,N)*DYF(J-1)+CTH(I,K,J-1,N)*DYF(J))/(2.d0*DY(J))
           END DO
         END DO
        END DO
       ELSEIF(N_TH == 2) THEN
!changes made for slope angles in the along the slope(X direction) and wall normal direction(Y direction)
!this part was running
       DO J=JSTART,JEND
         DO K=0,TNKZ
          DO I=0,NKX
             CF1(I,K,J)=CF1(I,K,J) 
     :                 + Ratio_gr_sc(N)*sin(ANG_BETA)*CTH(I,K,J,N) 
          END DO
         END DO
        END DO       

        DO J=2,NY
         DO K=0,TNKZ
           DO I=0,NKX
! Use second order interpolation
             CF2(I,K,J)=CF2(I,K,J) 
     :       - Ratio_gr_sc(N)*cos(ANG_BETA)
     :       *(CTH(I,K,J,N)*DYF(J-1)+CTH(I,K,J-1,N)*DYF(J))/(2.d0*DY(J))
           END DO
         END DO
        END DO    
       ENDIF

      ELSE
      DO J=2,NY
        DO K=1,TNKZ
          DO I=1,NKX
! Use second order interpolation
             CF2(I,K,J)=CF2(I,K,J)+RI_TAU(N)*
     &      (CTH(I,K,J,N)*DYF(J-1)+CTH(I,K,J-1,N)*DYF(J))/(2.d0*DY(J))
          END DO
        END DO
        K=0
        DO I=1,NKX 
             CF2(I,K,J)=CF2(I,K,J)+RI_TAU(N)*
     &      (CTH(I,K,J,N)*DYF(J-1)+CTH(I,K,J-1,N)*DYF(J))/(2.d0*DY(J))
        END DO
        I=0
        DO K=1,TNKZ
             CF2(I,K,J)=CF2(I,K,J)+RI_TAU(N)*
     &      (CTH(I,K,J,N)*DYF(J-1)+CTH(I,K,J-1,N)*DYF(J))/(2.d0*DY(J))
        END DO
      END DO 
      ENDIF

      ENDIF
! Now, compute the RHS vector for the scalar equations
! Since TH is defined at horizontal velocity points, the
! scalar update equation will be very similar to the horizontal
! velocity update.

! We will store the RHS scalar terms in CRTH, RTH
! The k-1 term for the R-K stepping is saved in FTH, CFTH

! First, build the RHS vector, use CRTH
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,TNKZ
          DO I=0,NKX
          CRTH(I,K,J,N)=CTH(I,K,J,N)
         ENDDO
       END DO
      END DO
! Add term from k-2 step to free up CFTH variable
      IF (RK_STEP .GT. 1) THEN
        DO J=JSTART_TH(N),JEND_TH(N)
          DO K=0,TNKZ
            DO I=0,NKX
              CRTH(I,K,J,N)=CRTH(I,K,J,N)+TEMP3*CFTH(I,K,J,N)
            END DO
          END DO
        END DO
       END IF

! Now that we have used the old CFTH, reset CFTH to zero
       DO J=0,NY+1
         DO K=0,TNKZ
           DO I=0,NKX
             CFTH(I,K,J,N)=0.d0
           END DO
         END DO
       END DO

! Now compute the explicit R-K term Ai
! Compile terms of Ai in CFi which will be saved for next time step
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,TNKZ
          DO I=0,NKX
            CFTH(I,K,J,N)=-(NU/PR(N)) * KX2(I) * CTH(I,K,J,N)
     &            - (NU/PR(N)) * KZ2(K) * CTH(I,K,J,N)
          END DO
        END DO
      END DO

C Add sponge in the theta field for open boundary condition
      IF (F_TYPE.EQ.4) THEN
       DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,TNKZ
          DO I=0,NKX
            CFTH(I,K,J,N)=CFTH(I,K,J,N) +(CU1(I,K,J)*sin(ANG_BETA) +
     :         0.5*(CU2(I,K,J+1)+CU2(I,K,J))*cos(ANG_BETA))*Grad_rho(N) ;
          END DO
        END DO
      END DO
      ELSE
      write(606,*) Grad_rho(N)
       IF(DEV_BACK_TH) THEN

        IF(N_TH == 1) THEN
         DO J=JSTART_TH(N),JEND_TH(N)
          DO K=0,TNKZ
           DO I=0,NKX
            CFTH(I,K,J,N)=CFTH(I,K,J,N) +
     :         0.5*(CU2(I,K,J+1)+CU2(I,K,J))*Grad_rho(N) ;
           END DO
          END DO
         END DO
        ELSEIF(N_TH == 2) THEN
!        write(6,*) 'Stratifiaction modification'changes of sign   
         DO J=JSTART_TH(N),JEND_TH(N)
          DO K=0,TNKZ
           DO I=0,NKX
            CFTH(I,K,J,N)=CFTH(I,K,J,N) + CU1(I,K,J)*sin(ANG_BETA)
     :         *Grad_rho(N) - 0.5*(CU2(I,K,J+1)+CU2(I,K,J))*
     :    cos(ANG_BETA)*Grad_rho(N);
           END DO
          END DO
         END DO
        ENDIF

       ENDIF
      ENDIF

      IF(SPONGE_FLAG)THEN
       IF ( GYF(NYC) .GE. y_sp_min ) THEN
         CALL sponge_th(N)
       ENDIF
      ENDIF
C End do number of passive scalars (N_TH)
      END DO

C   Add sponge in for open boundary condition

      IF(SPONGE_FLAG)THEN
       IF ( GYF(NYC) .GE. y_sp_min ) THEN
         CALL sponge_vel
       ENDIF
      ENDIF 

!     ADD FORCING TERM
      IF  (F_TYPE .eq. 4) THEN
        call forcing_vel  ! switching vel forcing 
        call forcing_th   ! switching fresh water cycle
      ENDIF 
  
C If we are considering an LES, then add the subgrid scale stress:
C Here, velocity and CFi should be in Fourier space
C The subgrid scale stress is added to CFi:   CFi=CFi - d/dx_i tau_ij

c      write(400,*) time_step, CU1(0,0,2), 10000
c      write(410,*) time_step, CU1(0,0,NY+1), 20000
      IF (MELTING_MODEL)THEN
       IF (RANK .eq. 0) THEN
        do N=1,N_TH
         dthdz_mean(n) = 0.0d0 ;
         dthdz_mean(n)=(dble(CTH(0,0,2,n))-dble(CTH(0,0,1,n)))
     &                    /(GYF(2)-GYF(1))
!                     (0.125*(CJOB_21(K,J,1)+CJOB_21(K,J+1,1) &
!                     + CJOB_21(K,J,2)+CJOB_21(K+1,J,2) )*    &
!                 dble(CTHX(0,k,j+1,n)-CTHX(0,k,j-1,n))       &
!               +  0.125*(CJOB_11(K,J,1)+CJOB_11(K,J+1,1)     &
!                     + CJOB_11(K,J,2)+CJOB_11(K+1,J,2) )*    &
!                  dble(CTHX(0,k+1,j,n)-CTHX(0,k-1,j,n)))/    &
!                  INT_JACOB(K,J)
         enddo
         write(666,166) time, dble(CU2(0,0,1)),dble(CTH(0,0,1,1)), 
     &          dble(CTH(0,0,1,2)),dthdz_mean(1), dthdz_mean(2),
     &          (C_sp_heat/L_heat)*(NU/Pr(1))*(1.0d0/0.92) 
     &          *dthdz_mean(1)
   
       ENDIF
       ENDIF

166    format (8f17.9)


      IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.20))) THEN
C If we have created new flow with random perturbations, wait for a
C spinup before applying the subgrid model for stability purposes
C In the process, Ui is converted to physical space
!          write(*,*) 'calling les here', rank

c        IF (USE_MPI) THEN
c           CALL  GHOST_ZERO_CHAN_MPI
c        ENDIF

          call les_chan
       
!           write(*,*) 'calling les done' ,rank
!          do j=0,NY+1
!          do k=0,NZM
!          do i=0,NXM
!            NU_T(I,K,J)=0.d0
!          end do
!          end do 
!          enddo
        
c          write(*,*) 'less is acting'
C Add the subgrid scale scalar flux to the scalar equations
          DO N=1,N_TH
            call les_chan_th(N)
          
!          do j=0,NY+1
!          do k=0,NZM
!          do i=0,NXM
!            KAPPA_T(I,K,J,N)=0.d0
!          end do
!          end do
!          end do

c          CALL FFT_XZ_TO_PHYSICAL(CTH(0,0,0,N),TH(0,0,0,N),0,NY+1)      
          END DO
C Now, we need to reapply the boundary conditions used for the momentum
C equations
          IF (USE_MPI) THEN
            CALL APPLY_BC_ANWM_MPI
          ELSE
            call APPLY_BC_ANWM_LOWER
            call APPLY_BC_ANWM_UPPER
          END IF 
 
C If we are using the wall model, make sure that the eddy viscosity is
C zero at the wall
         IF (USE_MPI) THEN
            IF ((U_BC_YMIN.eq.3).or.(W_BC_YMIN.eq.3)) THEN
            CALL APPLY_NUT_NWM_MPI_LOWER
            END IF 
            IF ((U_BC_YMAX.eq.3).or.(W_BC_YMAX.eq.3)) THEN
            CALL APPLY_NUT_NWM_MPI_UPPER
            END IF
         ELSE

           IF ((U_BC_YMIN.eq.3).or.(W_BC_YMIN.eq.3)) THEN
            DO I=0,NXM 
            DO K=0,NZM
              NU_T(I,K,2)=0.d0 
              DO N=1,N_TH
                KAPPA_T(I,K,2,N)=0.d0
              END DO
            END DO
            END DO
          END IF

           IF ((U_BC_YMAX.eq.3).or.(W_BC_YMAX.eq.3)) THEN
            DO I=0,NXM
            DO K=0,NZM
              NU_T(I,K,NY)=0.d0
              DO N=1,N_TH
                KAPPA_T(I,K,NY,N)=0.d0
              END DO
            END DO
            END DO
           END IF

        END IF


c       write(*,*) 'no less ' 

      ELSE 
C If the subgrid model hasn't been called, then it is necessary to 
C convert to physical space.
      CALL FFT_XZ_TO_PHYSICAL(CU1(0,0,-1),U1(0,0,-1),0,NY+1)
      CALL FFT_XZ_TO_PHYSICAL(CU2(0,0,-1),U2(0,0,-1),0,NY+1)
      CALL FFT_XZ_TO_PHYSICAL(CU3(0,0,-1),U3(0,0,-1),0,NY+1)

      CALL FFT_XZ_TO_PHYSICAL(CU1(0,0,NY+1),U1(0,0,NY+1),0,1)
      CALL FFT_XZ_TO_PHYSICAL(CU2(0,0,NY+1),U2(0,0,NY+1),0,1)
      CALL FFT_XZ_TO_PHYSICAL(CU3(0,0,NY+1),U3(0,0,NY+1),0,1)

!        CALL FFT_XZ_TO_PHYSICAL(CU1(0,0,0),U1(0,0,0),0,NY+1)
!        CALL FFT_XZ_TO_PHYSICAL(CU2(0,0,0),U2(0,0,0),0,NY+1)
!        CALL FFT_XZ_TO_PHYSICAL(CU3(0,0,0),U3(0,0,0),0,NY+1)
! Transform THETA to physical space for computation of nonlinear terms
! Here pass the first location in memory of the array for scalar n
        DO N=1,N_TH
          CALL FFT_XZ_TO_PHYSICAL(CTH(0,0,-1,N),TH(0,0,-1,N),0,NY+1)
          CALL FFT_XZ_TO_PHYSICAL(CTH(0,0,NY+1,N),TH(0,0,NY+1,N),0,1)
        END DO
      
      END IF

C Compute the nonlinear products in physical space, then transform
C back to Fourier space to compute the derivative.
C Here, we compute the horizontal derivatives of the nonlinear terms
C which will be treated with RKW3.  Nonlinear terms with vertical
C derivatives will be treated with Crank-Nicolson later
C Do terms one at a time to save on memory
C U1*U3
      DO J=JSTART,JEND
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,J)=U3(I,K,J)*U1(I,K,J)
          END DO
        END DO
      END DO
      
      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
      
      DO J=JSTART,JEND
        DO K=0,TNKZ
          DO I=0,NKX
            CF1(I,K,J)=CF1(I,K,J) - CIKZ(K) * CS1(I,K,J) 
            CF3(I,K,J)=CF3(I,K,J) - CIKX(I) * CS1(I,K,J) 
          END DO
        END DO
      END DO

C U1*U1
      DO J=JSTART,JEND
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,J)=U1(I,K,J)*U1(I,K,J)
          END DO
        END DO
      END DO
      
      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
      
      DO J=JSTART,JEND
        DO K=0,TNKZ
          DO I=0,NKX
            CF1(I,K,J)=CF1(I,K,J) - CIKX(I) * CS1(I,K,J) 
          END DO
        END DO
      END DO

C U3*U3
      DO J=JSTART,JEND
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,J)=U3(I,K,J)*U3(I,K,J)
          END DO
        END DO
      END DO
      
      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
      
      DO J=JSTART,JEND
        DO K=0,TNKZ
          DO I=0,NKX
            CF3(I,K,J)=CF3(I,K,J) - CIKZ(K) * CS1(I,K,J) 
          END DO
        END DO
      END DO


C U1*U2
      DO J=2,NY
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,J)=((DYF(J)*U1(I,K,J)
     &                +DYF(J-1)*U1(I,K,J-1))/(2.*DY(J))) 
     &                *U2(I,K,J)
          END DO
        END DO
      END DO
      
      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
      
      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NKX
            CF2(I,K,J)=CF2(I,K,J) - CIKX(I) * CS1(I,K,J) 
          END DO
        END DO
      END DO

C U3*U2
      DO J=2,NY
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,J)=((DYF(J)*U3(I,K,J)
     &                +DYF(J-1)*U3(I,K,J-1))/(2.*DY(J))) 
     &                *U2(I,K,J)
          END DO
        END DO
      END DO
      
      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
      
      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NKX
            CF2(I,K,J)=CF2(I,K,J) - CIKZ(K) * CS1(I,K,J)
          END DO
        END DO
      END DO

! Add the vertical derivative term explicitly
      DO J=JSTART,JEND
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,J)= ((U1(I,K,J+1)*U2(I,K,J+1) + U1(I,K,J)*U2(I,K,J+1)
     &     - U1(I,K,J)*U2(I,K,J) - U1(I,K,J-1)*U2(I,K,J))/(2.d0*DYF(J)))
          END DO
        END DO
      END DO

      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)

      DO J=JSTART,JEND
        DO K=0,TNKZ
          DO I=0,NKX
            CF1(I,K,J)=CF1(I,K,J) - CS1(I,K,J) 
         END DO
        END DO
      END DO
     
! Add the vertical derivative term explicitly
      DO J=JSTART,JEND
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,J)=((U3(I,K,J+1)*U2(I,K,J+1) + U3(I,K,J)*U2(I,K,J+1)
     &     - U3(I,K,J)*U2(I,K,J) - U3(I,K,J-1)*U2(I,K,J))/(2.d0*DYF(J)))
          END DO
        END DO
      END DO

      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)

      DO J=JSTART,JEND
        DO K=0,TNKZ
          DO I=0,NKX
            CF3(I,K,J)=CF3(I,K,J) - CS1(I,K,J) 
          END DO
        END DO
      END DO

! Add the vertical derivative term explicitly
      DO J=2,NY
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,J)=
     &    (0.25d0*(U2(I,K,J)+U2(I,K,J+1))**2.d0
     &    -0.25d0*(U2(I,K,J)+U2(I,K,J-1))**2.d0)/DY(J)
          END DO
        END DO
      END DO

      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)

      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NKX
            CF2(I,K,J)=CF2(I,K,J) - CS1(I,K,J)
          END DO
        END DO
      END DO


C -- At this point, we are done computing the nonlinear terms --

C Finally, Add CFi to CRi

      DO J=JSTART,JEND
        DO K=0,TNKZ
          DO I=0,NKX
            CR1(I,K,J)=CR1(I,K,J) + TEMP5 * CF1(I,K,J)
            CR3(I,K,J)=CR3(I,K,J) + TEMP5 * CF3(I,K,J)
          END DO
        END DO
      END DO
      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NKX
            CR2(I,K,J)=CR2(I,K,J) + TEMP5 * CF2(I,K,J)
          END DO
        END DO
      END DO

c      do j=JSTART,JEND
c        WRITE(300,*) j,CR1(0,0,J)
c      end do
!      write(*,*) gy(NY+1),gy(1) 
C Convert RHS terms to physical space
      CALL FFT_XZ_TO_PHYSICAL(CR1,R1,0,NY+1)                 
      CALL FFT_XZ_TO_PHYSICAL(CR2,R2,2,NY)                 
      CALL FFT_XZ_TO_PHYSICAL(CR3,R3,0,NY+1)                 

C Compute the vertical viscous term in physical space and add to RHS
C This is the explicit part of the Crank-Nicolson term
      DO J=JSTART,JEND
        DO K=0,NZM
          DO I=0,NXM
            R1(I,K,J)=R1(I,K,J)+TEMP1*
     &        (  ((U1(I,K,J+1) - U1(I,K,J)) / DY(J+1)  
     &           -(U1(I,K,J)   - U1(I,K,J-1)) / DY(J)) /DYF(J)  )
            R3(I,K,J)=R3(I,K,J)+TEMP1*
     &        (  ((U3(I,K,J+1) - U3(I,K,J)) / DY(J+1) 
     &           -(U3(I,K,J)   - U3(I,K,J-1)) / DY(J)) /DYF(J)  )
          END DO
        END DO
      END DO
      DO J=2,NY 
        DO K=0,NZM
          DO I=0,NXM
            R2(I,K,J)=R2(I,K,J)+TEMP1*
     &        (  ((U2(I,K,J+1) - U2(I,K,J))  / DYF(J) 
     &           -(U2(I,K,J)   - U2(I,K,J-1))/ DYF(J-1))/DY(J)  )
          END DO
        END DO
      END DO

C If we are using a subgrid model, add the eddy viscosity term
C This is an added viscosity that will be treated just like the 
C molecular viscosity with Crank-Nicolson for the vertical derivatives
c      IF (LES) then

      IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.40))) THEN  
C Note, NU_T is defined at GY points
      DO J=JSTART,JEND
        DO K=0,NZM
          DO I=0,NXM
            R1(I,K,J)=R1(I,K,J)+TEMP2*
     &        (  (NU_T(I,K,J+1) * (U1(I,K,J+1) - U1(I,K,J)) / DY(J+1)  
     &         -  NU_T(I,K,J) * (U1(I,K,J)   - U1(I,K,J-1)) / DY(J))
     &               /DYF(J)  )
            R3(I,K,J)=R3(I,K,J)+TEMP2*
     &        (  (NU_T(I,K,J+1) * (U3(I,K,J+1) - U3(I,K,J)) / DY(J+1) 
     &        - NU_T(I,K,J) * (U3(I,K,J)   - U3(I,K,J-1)) / DY(J)) 
     &              /DYF(J)  )
          END DO
        END DO
      END DO
! Here, interpolate NU_T to GYF points
      DO J=2,NY 
        DO K=0,NZM
          DO I=0,NXM
            R2(I,K,J)=R2(I,K,J)+TEMP2_v*
     &     ((0.5d0*(NU_T(I,K,J)+NU_T(I,K,J+1))*(U2(I,K,J+1)-U2(I,K,J))
     &                                              / DYF(J) 
     &    -0.5d0*(NU_T(I,K,J)+NU_T(I,K,J-1))*(U2(I,K,J)-U2(I,K,J-1))
     &                                          / DYF(J-1))   /DY(J)  )
          END DO
        END DO
      END DO
      END IF


C -- Here, we are done with computation of Velocity RHS, explicit terms --

C Now, build the explicit RHS terms for the passive scalar(s)

      DO N=1,N_TH
! Do for each scalar:

! Compute the nonlinear terms that are present in the explicit term A
! U1*TH
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,J)=TH(I,K,J,N)*U1(I,K,J)
          END DO
        END DO
      END DO
      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,TNKZ
          DO I=0,NKX
            CFTH(I,K,J,N)=CFTH(I,K,J,N) - CIKX(I) * CS1(I,K,J)
          END DO
        END DO
      END DO
! U3*TH 
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,J)=TH(I,K,J,N)*U3(I,K,J)
          END DO
        END DO
      END DO
      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,TNKZ
          DO I=0,NKX
            CFTH(I,K,J,N)=CFTH(I,K,J,N) - CIKZ(K) * CS1(I,K,J)
          END DO
        END DO
      END DO
! U2*TH
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,NZM
          DO I=0,NXM
         S1(I,K,J)=((TH(I,K,J+1,N)*U2(I,K,J+1) + TH(I,K,J,N)*U2(I,K,J+1)
     &    -TH(I,K,J,N)*U2(I,K,J)-TH(I,K,J-1,N)*U2(I,K,J))/(2.d0*DYF(J)))
          END DO
        END DO
      END DO

c      write(400,*)time_step,sum(U2(0:NXM,0:NZM,1)), 1000, JSTART_TH(N)
c      write(410,*)time_step,sum(U2(0:NXM,0:NZM,NY)),2000,JEND_TH(N)

      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,TNKZ
          DO I=0,NKX
            CFTH(I,K,J,N)=CFTH(I,K,J,N) - CS1(I,K,J)
          END DO
        END DO
      END DO

! We are done with the horizontal derivatives of the nonlinear terms

! Add CFTH to the RHS vector CRTH
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,TNKZ
          DO I=0,NKX
            CRTH(I,K,J,N)=CRTH(I,K,J,N) + TEMP5 * CFTH(I,K,J,N)
          END DO
        END DO
      END DO
! Done with computation of RHS, explicit terms for the THETA equation
! Transform back to physical space

      CALL FFT_XZ_TO_PHYSICAL(CRTH(0,0,0,N),RTH(0,0,0,N),0,NY+1)    

! Compute the Explicit part of the Crank-Nicolson terms for the TH equation
! First, the vertical derivative viscous term
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,NZM
          DO I=0,NXM
            RTH(I,K,J,N)=RTH(I,K,J,N)+(TEMP1/PR(N))*(
     &            ((TH(I,K,J+1,N) - TH(I,K,J,N)) / DY(J+1)
     &            -(TH(I,K,J,N) - TH(I,K,J-1,N)) / DY(J)) / DYF(J) )
          END DO
        END DO
      END DO

c      write(400,*)time_step,DY(1), TEMP1,PR(N)
c      write(410,*)time_step,DY(NY),TEMP1,PR(N)

! If we are using a subgrid model (LES) then add the eddy diffusivity here
! Note, KAPPA_T is defined at GY points
c      IF (LES) THEN

      IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.40))) THEN 
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,NZM
          DO I=0,NXM
            RTH(I,K,J,N)=RTH(I,K,J,N)+TEMP2*(
     &     (KAPPA_T(I,K,J+1,N)*(TH(I,K,J+1,N)-TH(I,K,J,N))/DY(J+1)
     &     -KAPPA_T(I,K,J,N)*(TH(I,K,J,N)-TH(I,K,J-1,N))/DY(J))/DYF(J))
          END DO
        END DO
      END DO  
      END IF



C -- Now, timestep the passive scalar equation --
C      We solve the the passive scalar before the velocity so that
C      it is advected with the velocity from the previous R-K step
C      which we have already made divergence free 
 
! Solve the implicit equation for THETA
! Note that the system size is NY+1, but only 1..NY are used

! Initialize the matrix used to store implicit coefficients
      DO J=0,NY+1
        DO I=0,NXM
          MATL(I,J)=0.
          MATD(I,J)=1.
          MATU(I,J)=0.
          VEC(I,J)=0.
        END DO
      END DO 
      
 
c      do j=0,NY+1
c      do k=0,NZM
c      do i=0,NXM
c        th_r(i,k,j) = U2(i,k,j)
c      end do
c      end do
c      enddo
c      
c      sum_th=0.d0
c      j=NY
c      do k=0,NZM
c      do i=0,NXM
c        sum_th = sum_th + th_R(i,k,j)
c      end do
c      end do
c      thrms_ny=sum_th/(float(NZ)*float(NX))
c
c      sum_th=0.d0
c      j=1
c      do k=0,NZM
c      do i=0,NXM
c        sum_th = sum_th + th_R(i,k,j)
c      end do
c      end do
c      thrms_1=sum_th/(float(NZ)*float(NX))

c      CALL FFT_XZ_TO_FOURIER(TH(0,0,-1,n),CTH(0,0,-1,n),0,NY+1)
c      CALL FFT_XZ_TO_FOURIER(TH(0,0,NY+1,n),CTH(0,0,NY+1,n),0,1) 
c       call GHOST_CHAN_MPI
c       call GHOST_ZERO_CHAN_MPI
c      sum_th=0.0
c      j=NY
c     do k=0,NZM
c      do i=0,NXM
c        sum_th = sum_th + (abs(th_R(i,k,j)-dble(CU2(0,0,j))))**2.
c      end do
c      end do
c      thrms_ny=sqrt(sum_th/(float(NZ)*float(NX)))
c      j=1
c      sum_th=0.0
c      do k=0,NZM
c      do i=0,NXM
c        sum_th=sum_th+(abs(th_R(i,k,j)-dble(CU2(0,0,j))))**2.
c      end do
c      end do
      
c       thrms_1=sqrt(sum_th/(float(NZ)*float(NX))) 
       
c      write(400,*) time_step, rk_step, CTH(5,5,1,1),1000
c      write(410,*) time_step, rk_step, CTH(5,5,NY,1),2000
c
c        write(400,*) time_step, rk_step, CTH(0,0,-1,n),1000
c        write(410,*) time_step, rk_step, CTH(0,0,NY-2,n),1000

c          CALL FFT_XZ_TO_PHYSICAL(CTH(0,0,-1,n),TH(0,0,-1,n),0,NY+1)
c          CALL FFT_XZ_TO_PHYSICAL(CTH(0,0,NY+1,n),TH(0,0,NY+1,n),0,1)
 
c       sum_th=0.0
c      j=NY+1
c      do k=0,NZM
c      do i=0,NXM
c        sum_th = sum_th + TH(i,k,j,n)
c      end do
c      end do
c      thrms_ny=sum_th/(float(NZ)*float(NX))
c      sum_th=0.0
c      j=2
c      do k=0,NZM
c      do i=0,NXM
c        sum_th = sum_th + TH(i,k,j,n)
c      end do
c      end do
c      thrms_1=sum_th/(float(NZ)*float(NX))
c
c      write(400,*) time_step, rk_step, thrms_1,1000
c      write(410,*) time_step, rk_step, thrms_ny,2000

    
! Build implicit matrix
! Use quasi-second order interpolation for TH on GY points
      DO K=0,NZM
        DO J=JSTART_TH(N),JEND_TH(N)
          DO I=0,NXM
            MATL(I,J) = -(TEMP1/PR(N)) / (DY(J)*DYF(J))
            MATD(I,J) = 1. + (TEMP1/PR(N)) / (DY(J+1)*DYF(J))
     &           +(TEMP1/PR(N)) / (DY(J)*DYF(J))
            MATU(I,J)=-(TEMP1/PR(N)) / (DY(J+1)*DYF(J))
! Define RHS vector
            VEC(I,J)=RTH(I,K,J,N)
          END DO
        END DO
! IF using a subgrid model (LES) then add the eddy diffusivity part implicitly
c        IF (LES) THEN
      IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.40))) THEN
c        DO J=JSTART_TH(N),JEND_TH(N)
        DO J=0,NY-1

          DO I=0,NXM   
            MATL(I,J) = MATL(I,J) - TEMP2 * KAPPA_T(I,K,J,N) 
     &                                 / (DY(J)*DYF(J))
            MATD(I,J) = MATD(I,J)+ TEMP2 * KAPPA_T(I,K,J+1,N)
     &                                 / (DY(J+1)*DYF(J))
     &                           + TEMP2 * KAPPA_T(I,K,J,N)
     &                                 / (DY(J)*DYF(J))
            MATU(I,J) = MATU(I,J)- TEMP2 * KAPPA_T(I,K,J+1,N)
     &                                / (DY(J+1)*DYF(J))
          END DO
        END DO
        END IF

! If we are using MPI, then solve the implicit system in separate forward
! and backward sweeps for efficiency
          IF (USE_MPI) THEN
             CALL APPLY_BC_TH_MPI(MATL,MATD,MATU,VEC,K,N)
! If we are using MPI, split the implicit solve into foward and
! backward sweeps for efficiency
             CALL THOMAS_FORWARD_REAL_MPI(MATL,MATD,MATU,VEC,NY,NX)
             CALL THOMAS_BACKWARD_REAL_MPI(MATL,MATD,MATU,VEC,NY,NX)
          ELSE
! Else we are running in serial mode
             CALL APPLY_BC_TH_LOWER(MATL,MATD,MATU,VEC,N)
             CALL APPLY_BC_TH_UPPER(MATL,MATD,MATU,VEC,N)
             CALL THOMAS_REAL(MATL,MATD,MATU,VEC,NY+1,NXM)
! Apply the boundary conditions to our linear system
          END IF

        DO J=JSTART_TH(N),JEND_TH(N)
          DO I=0,NXM
            TH(I,K,J,N)=VEC(I,J)
          END DO
        END DO

! END do k
      END DO
     
!      if(rank == 0) then
!      write(6,*)'TH down = ', TH(0,0,2,N)
!      end if 
!      if(rank == 3) then
!      write(6,*)'TH up = ', TH(0,0,NY-1,N)
!      end if
c      sum_th=0.0
c      j=NY+1
c      do k=0,NZM
c      do i=0,NXM
c        sum_th = sum_th + TH(i,k,j,n)
c      end do
c      end do
c      thrms_ny=sum_th/(float(NZ)*float(NX)) 

c      sum_th=0.0
c      j=NY-1
c      do k=0,NZM
c      do i=0,NXM
c        sum_th = sum_th + TH(i,k,j,n)
c      end do
c      end do
c      thrms_1=sum_th/(float(NZ)*float(NX))
   

c      write(410,*)JEND_TH(N)+1,GYF(NY+1),thrms_ny,thrms_1,
c     &(thrms_ny-thrms_1)/(2.*DYF(NY)),(GYF(NY+1)-GYF(NY-1))/(2.*DYF(NY))
       


c        write(400,*)  TIME_STEP, sum(TH(1:NX,1:NZ,1,N)), 1000
c        write(410,*)  TIME_STEP, sum(TH(1:NX,1:NZ,NY,N)), 2000

! End do number of passive scalars
        END DO

        IF (RANK == 0) THEN
!          write(6,*) 'Here', JSTART_TH(1), JSTART_TH(2), 
!     &    DY(1)*TH_BC_YMIN_C1(2)*
!     &               (5.0d0*PR(2)/PR(1))*(C_sp_heat/L_heat)*
!     &               (TH(1,1,1,1) + TH_BACK(1,1))*dthdz_mean(1)/m_FP


          DO I=0,NXM
           WRITE(656,*) GX(I), TH(I,0,1,1), TH(I,0,1,2) 
          ENDDO
          close(656)
        ENDIF         
        

c      FIRST_TIME=.FALSE.
! Save statistics to an output file
c      IF (MOD(TIME_STEP,SAVE_STATS_INT).EQ.0) THEN
c       write(*,*) 'hello i am fool'
C Convert velocity back to Fourier space
c      call fft_xz_to_fourier(U1(0,0,0),CU1(0,0,0),0,NY+1)
c      call fft_xz_to_fourier(U2(0,0,0),CU2(0,0,0),0,NY+1)
c      call fft_xz_to_fourier(U3(0,0,0),CU3(0,0,0),0,NY+1)
c      DO N=1,N_TH
c       CALL FFT_XZ_TO_FOURIER(TH(0,0,0,N),CTH(0,0,0,N),0,NY+1)
c      END DO
c            CALL SAVE_STATS(.FALSE.)
C convert to physical space.
c        CALL FFT_XZ_TO_PHYSICAL(CU1(0,0,0),U1(0,0,0),0,NY+1)
c        CALL FFT_XZ_TO_PHYSICAL(CU2(0,0,0),U2(0,0,0),0,NY+1)
c        CALL FFT_XZ_TO_PHYSICAL(CU3(0,0,0),U3(0,0,0),0,NY+1)
! Transform THETA to physical space for computation of nonlinear terms
! Here pass the first location in memory of the array for scalar n
c        DO N=1,N_TH
c          CALL FFT_XZ_TO_PHYSICAL(CTH(0,0,0,N),TH(0,0,0,N),0,NY+1)
c        END DO
c      ENDIF



C If a wall model is being used, compute the value of the wall gradient
C This will be needed in the next step when solving the implicit systems
C for the horizontal velocity components
      IF (USE_MPI) THEN
         CALL APPLY_WALL_MODEL_MPI
      ELSE
        IF ((U_BC_YMIN.EQ.3).or.(W_BC_YMIN.eq.3)) THEN
         CALL WALL_MODEL_LOWER
        END IF
        IF ((U_BC_YMAX.EQ.3).or.(W_BC_YMAX.eq.3)) THEN
         CALL WALL_MODEL_UPPER
        END IF
      ENDIF



C Initialize the matrix to zeros to be used for implicit solves
C Note that the system size is NY+1, but only 1..NY are used

! Initialize the matrix used to store implicit coefficients
      DO J=0,NY+1
        DO I=0,NXM
          MATL(I,J)=0.
          MATD(I,J)=1.
          MATU(I,J)=0.
          VEC(I,J)=0.
        END DO
      END DO 

C Build implicit matrix for U2
      DO K=0,NZM
        DO J=2,NY
          DO I=0,NXM
            MATL(I,J)= -TEMP1/(DYF(J-1)*DY(J))
            MATD(I,J)=1.+TEMP1/(DYF(J)*DY(J)) + TEMP1/(DYF(J-1)*DY(J)) 
            MATU(I,J)= -TEMP1/(DYF(J)*DY(J))
            VEC(I,J)=R2(I,K,J)
          END DO 
        END DO
c        IF (LES) THEN
       IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.40))) THEN
! IF using a subgrid model (LES) then add the eddy viscosity part implicitly
        DO J=2,NY
          DO I=0,NXM
            MATL(I,J) = MATL(I,J) 
     &    -TEMP2_v* 0.5d0*(NU_T(I,K,J)+NU_T(I,K,J-1))/(DYF(J-1)*DY(J))
            MATD(I,J) = MATD(I,J) 
     &    +TEMP2_v* 0.5d0*(NU_T(I,K,J)+NU_T(I,K,J+1))/(DYF(J)*DY(J))
     &    + TEMP2_v* 0.5d0*(NU_T(I,K,J)+NU_T(I,K,J-1))/(DYF(J-1)*DY(J))
            MATU(I,J) = MATU(I,J) 
     &    - TEMP2_v* 0.5d0*(NU_T(I,K,J)+NU_T(I,K,J+1))/(DYF(J)*DY(J))
          END DO
        END DO
        END IF

        IF (USE_MPI) THEN
! First, apply the boundary conditions
          CALL APPLY_BC_U2_MPI(MATL,MATD,MATU,VEC,K)
         !    do j = 0,NY+1
         !     write(400,*) j,VEC(2,j)
         !    enddo 
! If we are using MPI, split the implicit solve into forward and
! backward sweeps for efficiency
          CALL THOMAS_FORWARD_REAL_MPI(MATL,MATD,MATU,VEC,NY,NX)
          CALL THOMAS_BACKWARD_REAL_MPI(MATL,MATD,MATU,VEC,NY,NX)
             
        ELSE
C Else, we are running in serial mode
C Set the boundary conditions for U2
          CALL APPLY_BC_2_LOWER(MATL,MATD,MATU,VEC,K)
          CALL APPLY_BC_2_UPPER(MATL,MATD,MATU,VEC,K)
C Now, solve the tridiagonal system for U2(i,:,k)
          CALL THOMAS_REAL(MATL,MATD,MATU,VEC,NY+1,NXM)
        END IF

        DO J=1,NY+1
          DO I=0,NXM
            U2(I,K,J)=VEC(I,J)
          END DO

        END DO
! End do k
      END DO 

c      write(400,*)  TIME_STEP, sum(U2(1:NX,1:NZ,1)), 1000
c      write(410,*)  TIME_STEP, sum(U2(1:NX,1:NZ,NY)), 2000 
           
C Solve for U1
C Note, here the matrix will be indexed from 1...NY+1 corresponding to U1(0:NY)

! Initialize the matrix used to store implicit coefficients
      DO J=0,NY+1
        DO I=0,NXM
          MATL(I,J)=0.
          MATD(I,J)=1.
          MATU(I,J)=0.
          VEC(I,J)=0.
        END DO
      END DO 

C Build the implicit system of equations for U1 
      DO K=0,NZM
        DO J=JSTART,JEND
          DO I=0,NXM
            MATL(I,J)=-TEMP1/(DY(J)*DYF(J))
            MATD(I,J)=1.-TEMP1*(-1./(DY(J+1)*DYF(J))
     &         -1./(DY(J)*DYF(J))) 
            MATU(I,J)=-TEMP1/(DY(J+1)*DYF(J))
            VEC(I,J)=R1(I,K,J)
          END DO
        END DO
! IF using a subgrid model (LES) then add the eddy viscosity part implicitly
c        IF (LES) THEN
       IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.40))) THEN
        DO J=JSTART,JEND
          DO I=0,NXM
            MATL(I,J) = MATL(I,J) - TEMP2 * NU_T(I,K,J) 
     &                               / (DY(J)*DYF(J))
            MATD(I,J) = MATD(I,J) + TEMP2 * NU_T(I,K,J+1)
     &                              / (DY(J+1)*DYF(J))
     &                            + TEMP2 * NU_T(I,K,J)
     &                              / (DY(J)*DYF(J))
            MATU(I,J) = MATU(I,J) - TEMP2 * NU_T(I,K,J+1)
     &                             / (DY(J+1)*DYF(J))
          END DO
        END DO
        END IF

        IF (USE_MPI) THEN
! First, apply the boundary conditions
          CALL APPLY_BC_U1_MPI(MATL,MATD,MATU,VEC,K)
! If we are using MPI, split the implicit solve into forward and
! backward sweeps for efficiency
          CALL THOMAS_FORWARD_REAL_MPI(MATL,MATD,MATU,VEC,NY,NX)
          CALL THOMAS_BACKWARD_REAL_MPI(MATL,MATD,MATU,VEC,NY,NX)
        ELSE
C Else, we are running in serial mode
C Set the boundary conditions for U1
          CALL APPLY_BC_1_LOWER(MATL,MATD,MATU,VEC,K)
          CALL APPLY_BC_1_UPPER(MATL,MATD,MATU,VEC,K)
C Now, solve the tridiagonal system for U1(:,k,:)
          CALL THOMAS_REAL(MATL,MATD,MATU,VEC,NY+1,NXM)
        END IF

        DO J=JSTART-1,JEND+1
          DO I=0,NXM
            U1(I,K,J)=VEC(I,J)
          END DO
        END DO
! End do k
      END DO

c      write(400,*)  TIME_STEP, sum(U1(1:NX,1:NZ,2)), 1000
c      write(410,*)  TIME_STEP, sum(U1(1:NX,1:NZ,NY+1)), 2000 

! Initialize the matrix used to store implicit coefficients
      DO J=0,NY+1
        DO I=0,NXM
          MATL(I,J)=0.
          MATD(I,J)=1.
          MATU(I,J)=0.
          VEC(I,J)=0.
        END DO
      END DO 

C Solve for U3
C Note, here the matrix will be indexed from 1...NY+1 corresponding to U1(0:NY)
C Build the implicit system of equations for U3
      DO K=0,NZM
        DO J=JSTART,JEND
          DO I=0,NXM
            MATL(I,J)=-TEMP1/(DY(J)*DYF(J))
            MATD(I,J)=1.-TEMP1*(-1./(DY(J+1)*DYF(J))
     &         -1./(DY(J)*DYF(J)))
            MATU(I,J)=-TEMP1/(DY(J+1)*DYF(J))
            VEC(I,J)=R3(I,K,J)
          END DO
        END DO
! IF using a subgrid model (LES) then add the eddy viscosity part implicitly
c        IF (LES) THEN
      IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.40))) THEN
        DO J=JSTART,JEND
          DO I=0,NXM
            MATL(I,J) = MATL(I,J) - TEMP2 * NU_T(I,K,J)
     &                               / (DY(J)*DYF(J))
            MATD(I,J) = MATD(I,J) + TEMP2 * NU_T(I,K,J+1)
     &                              / (DY(J+1)*DYF(J))
     &                            + TEMP2 * NU_T(I,K,J)
     &                              / (DY(J)*DYF(J))
            MATU(I,J) = MATU(I,J) - TEMP2 * NU_T(I,K,J+1)
     &                             / (DY(J+1)*DYF(J))
          END DO
        END DO
        END IF

        IF (USE_MPI) THEN
! First, apply the boundary conditions
          CALL APPLY_BC_U3_MPI(MATL,MATD,MATU,VEC,K)
! If we are using MPI, split the implicit solve into forward and
! backward sweeps for efficiency
          CALL THOMAS_FORWARD_REAL_MPI(MATL,MATD,MATU,VEC,NY,NX)
          CALL THOMAS_BACKWARD_REAL_MPI(MATL,MATD,MATU,VEC,NY,NX)
        ELSE
C Else, we are running in serial mode
C Set the boundary conditions for U3
          CALL APPLY_BC_3_LOWER(MATL,MATD,MATU,VEC,K)
          CALL APPLY_BC_3_UPPER(MATL,MATD,MATU,VEC,K)
C Now, solve the tridiagonal system for U3(i,:,k)
          CALL THOMAS_REAL(MATL,MATD,MATU,VEC,NY+1,NXM)
        END IF

        DO J=JSTART-1,JEND+1
          DO I=0,NXM
            U3(I,K,J)=VEC(I,J)
          END DO
        END DO
! End do k
      END DO


c      FIRST_TIME=.FALSE.
! Save statistics to an output file
c      IF (MOD(TIME_STEP,SAVE_STATS_INT).EQ.0) THEN
c       write(*,*) 'hello i am fool'
C Convert velocity back to Fourier space
c      call fft_xz_to_fourier(U1(0,0,0),CU1(0,0,0),0,NY+1)
c      call fft_xz_to_fourier(U2(0,0,0),CU2(0,0,0),0,NY+1)
c      call fft_xz_to_fourier(U3(0,0,0),CU3(0,0,0),0,NY+1)
c      DO N=1,N_TH
c       CALL FFT_XZ_TO_FOURIER(TH(0,0,0,N),CTH(0,0,0,N),0,NY+1)
c      END DO
c            CALL SAVE_STATS(.FALSE.)
C convert to physical space.
c        CALL FFT_XZ_TO_PHYSICAL(CU1(0,0,0),U1(0,0,0),0,NY+1)
c        CALL FFT_XZ_TO_PHYSICAL(CU2(0,0,0),U2(0,0,0),0,NY+1)
c        CALL FFT_XZ_TO_PHYSICAL(CU3(0,0,0),U3(0,0,0),0,NY+1)
! Transform THETA to physical space for computation of nonlinear terms
! Here pass the first location in memory of the array for scalar n
c        DO N=1,N_TH
c          CALL FFT_XZ_TO_PHYSICAL(CTH(0,0,0,N),TH(0,0,0,N),0,NY+1)
c        END DO
c      ENDIF

C -- Done getting U1hat, U2hat, U3hat at new RK Step --

C If Variable timestepping and done with one full R-K step, update
C DELTA_T based on the specified CFL number
C This is not parallelized and should be used only in the serial
C version to ensure that each process uses the same timestep
C      IF ((.NOT.USE_MPI).and.(VARIABLE_DT).and.(RK_STEP.eq.3)
C     &        .and.(MOD(TIME_STEP,UPDATE_DT).EQ.0)) THEN
      IF ((VARIABLE_DT).and.(RK_STEP.eq.3)
     &        .and.(MOD(TIME_STEP,UPDATE_DT).EQ.0)) THEN
        CALL COURANT
      END IF
 
! Transform TH and U to Fourier Space 
      CALL FFT_XZ_TO_FOURIER(U1(0,0,0),CU1(0,0,0),0,NY+1)
      CALL FFT_XZ_TO_FOURIER(U2(0,0,0),CU2(0,0,0),0,NY+1)
      CALL FFT_XZ_TO_FOURIER(U3(0,0,0),CU3(0,0,0),0,NY+1)
      DO N=1,N_TH
        CALL FFT_XZ_TO_FOURIER(TH(0,0,0,N),CTH(0,0,0,N),0,NY+1)
      END DO
        
       
C Begin second step of the Fractional Step algorithm, making u divergence free
C The following subroutine projects Uhat onto divergence free space
       CALL GHOST_CHAN_MPI
       !write(400,*) time_step,sum(U1(0:NX,0:NZ,NY)), NKX,NX/3
       !write(410,*) time_step,sum(U1(0:NX,0:NZ,NY+1)) 
      CALL REM_DIV_CHAN
    
      !CALL FFT_XZ_TO_PHYSICAL(CU1,U1,0,NY+1)
        !write(430,*) time_step,sum(U1(0:NX,0:NZ,NY+1))
        !  CALL FFT_XZ_TO_FOURIER(U1,CU1,0,NY+1)
        !  CALL FFT_XZ_TO_PHYSICAL(CU1,U1,0,NY+1)
        !write(440,*) time_step,sum(U1(0:NX,0:NZ,NY+1))
      !CALL FFT_XZ_TO_FOURIER(U1,CU1,0,NY+1)
        
      
C Now, phi is stored in CR1, use this to update the pressure field
C Note, here we divide by H_BAR since it was absorbed into PHI in REM_DIV
      DO J=JSTART,JEND
        DO K=0,TNKZ
          DO I=0,NKX
            CP(I,K,J)=CP(I,K,J)+CR1(I,K,J)/TEMP4
          END DO
        END DO
      END DO

!      IF (MOD(TIME_STEP,1000).EQ.0) THEN
!        CALL POISSON_P_CHAN
!      ENDIF

      RETURN
      END
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE REM_DIV_CHAN
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      
C Compute varphi, store in variable CR1.
C Solves for phi in computational space
C H_BAR has been absorbed into PHI, so we are solving for H_BAR*PHI

      INCLUDE 'header'
      INTEGER I,J,K
 
C First, Initialize the matrix components
      DO J=0,NY+1
        DO I=0,NKX
          MATL_C(I,J)=0.
          MATD_C(I,J)=1.
          MATU_C(I,J)=0.
          VEC_C(I,J)=(0.,0.)
        END DO
      END DO

C The 2d FFT of Ui should have been taken and stored in CUi
C Solving for phi amounts to solving a tridiagonal system
C First, construct the system to be solved
      DO K=0,TNKZ
        DO J=1,NY
          DO I=0,NKX
            MATL_C(I,J)=1./(DY(J)*DYF(J))
            MATD_C(I,J)=-KX2(I)-KZ2(K)
     &         -1./(DY(J+1)*DYF(J))-1./(DY(J)*DYF(J))
            MATU_C(I,J)=1./(DY(J+1)*DYF(J))
          END DO
        END DO

C Now, create the RHS vector
        DO J=1,NY         
          DO I=0,NKX
            VEC_C(I,J)=(CIKX(I)*CU1(I,K,J) 
     &            + (CU2(I,K,J+1)-CU2(I,K,J))/DYF(J) 
     &            + CIKZ(K)*CU3(I,K,J))
          END DO
        END DO

        IF (USE_MPI) THEN
C If we are using the MPI package...
          CALL APPLY_BC_REM_DIV_MPI(MATL_C,MATD_C,MATU_C,VEC_C,K)
C First, do all forward sweeps
          CALL THOMAS_FORWARD_COMPLEX_MPI(MATL_C,MATD_C,MATU_C,VEC_C
     &                                 ,NY,NX/3)
C Now, do the backward sweeps
          CALL THOMAS_BACKWARD_COMPLEX_MPI(MATL_C,MATD_C,MATU_C,VEC_C
     &                                 ,NY,NX/3)
        ELSE
C Else we are running in serial mode
        DO I=0,NKX
          IF ((K.EQ.0).AND.(I.EQ.0)) THEN
C Use homogeneous dirichlet BCS for kx=kz=0 component at bottom wall
C Otherwise the matrix will be singular
            MATL_C(I,1)=0. 
            MATD_C(I,1)=1.
            MATU_C(I,1)=0.
            VEC_C(I,1)=(0.,0.)

            MATL_C(I,NY)=1.
            MATD_C(I,NY)=-1.
            MATU_C(I,NY)=0.
            VEC_C(I,NY)=(0.,0.)
          ELSE
C Use Dirichlet boundary conditions, dp/dz=0 at walls
            MATL_C(I,1)=0.
            MATD_C(I,1)=1.
            MATU_C(I,1)=-1.
            VEC_C(I,1)=(0.,0.)

            MATL_C(I,NY)=1.
            MATD_C(I,NY)=-1.
            MATU_C(I,NY)=0.
            VEC_C(I,NY)=(0.,0.)
          END IF
        END DO
C Now solve the tridiagonal system for phi, store in CR1
        CALL THOMAS_COMPLEX(MATL_C,MATD_C,MATU_C,VEC_C,NY,NKX)
        END IF


        DO J=1,NY
          DO I=0,NKX
            CR1(I,K,J)=VEC_C(I,J)
          END DO
        END DO

      END DO

! NEED TO EXCAHNGE VALUE OF CR1 BETWEEN PROCESSORS

       CALL GHOST_PHI_MPI
       
C Now, Solve for CUi, the divergenceless velocity field
      DO J=1,NY
C      DO J=JSTART,JEND
        DO K=0,TNKZ
          DO I=0,NKX
            CU1(I,K,J)=CU1(I,K,J)-CIKX(I)*CR1(I,K,J)
            CU3(I,K,J)=CU3(I,K,J)-CIKZ(K)*CR1(I,K,J)           
          END DO
        END DO
      END DO
      
      

      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NKX
            CU2(I,K,J)=CU2(I,K,J)-(CR1(I,K,J)
     &             -CR1(I,K,J-1))/DY(J)
          END DO
        END DO
      END DO

c      write(400,*) time_step, CU2(0,0,1)
c      write(410,*) time_step, CU2(0,0,NY)   

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE POISSON_P_CHAN
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C We have CUi, need to compute CP.  Solve tridiagonal system exactly

      INCLUDE 'header'
      INCLUDE "mpif.h"
      include 'header_mpi'
  
      INTEGER I,J,K,N
    
      if (rank .eq. 0) then
      WRITE(*,*) 'COMPUTING CP FROM CUI'
      endif
 
C First, construct the RHS vector, (dui/dxj)(duj/dxi) 
      DO J=2,NYM
        DO K=0,TNKZ
          DO I=0,NKX
            CF1(I,K,J)=CIKX(I)*CU1(I,K,J)
            CF2(I,K,J)=(CU2(I,K,J+1)-CU2(I,K,J))/DYF(J)
            CF3(I,K,J)=CIKZ(K)*CU3(I,K,J)
          END DO
        END DO
      END DO

      CALL FFT_XZ_TO_PHYSICAL(CF1,F1,0,NY+1)
      CALL FFT_XZ_TO_PHYSICAL(CF2,F2,0,NY+1)
      CALL FFT_XZ_TO_PHYSICAL(CF3,F3,0,NY+1)
      
      DO J=2,NYM
        DO K=0,NZM
          DO I=0,NXM
            F1(I,K,J)=F1(I,K,J)**2.
            F2(I,K,J)=F2(I,K,J)**2.
            F3(I,K,J)=F3(I,K,J)**2.
          END DO
        END DO
      END DO

      CALL FFT_XZ_TO_FOURIER(F1,CF1,0,NY+1)
      CALL FFT_XZ_TO_FOURIER(F2,CF2,0,NY+1)
      CALL FFT_XZ_TO_FOURIER(F3,CF3,0,NY+1)
      
C Now we have the diagonal terms, add to the rhs term
      DO J=2,NYM
        DO K=0,TNKZ
          DO I=0,NKX
            CS1(I,K,J)=CF1(I,K,J)+CF2(I,K,J)+CF3(I,K,J)
          END DO
        END DO
      END DO

C Now get the first of the off-diagonal terms
      DO J=2,NYM
        DO K=0,TNKZ
          DO I=0,NKX
            CF1(I,K,J)=(CU1(I,K,J+1)-CU1(I,K,J-1))/(2.*DYF(J))
            CF2(I,K,J)=CIKX(I)*0.5*(CU2(I,K,J)+CU2(I,K,J+1))
          END DO
        END DO
      END DO

      CALL FFT_XZ_TO_PHYSICAL(CF1,F1,0,NY+1)
      CALL FFT_XZ_TO_PHYSICAL(CF2,F2,0,NY+1)

C Compute product
      DO J=2,NYM
        DO K=0,NZM
          DO I=0,NXM
            F1(I,K,J)=2.*F1(I,K,J)*F2(I,K,J)
          END DO
        END DO
      END DO
      
      CALL FFT_XZ_TO_FOURIER(F1,CF1,0,NY+1)

C Add to RHS term
      DO J=2,NYM
        DO K=0,TNKZ
          DO I=0,NKX 
            CS1(I,K,J)=CS1(I,K,J)+CF1(I,K,J)
          END DO
        END DO
      END DO

C Now get the second of the off-diagonal terms
      DO J=2,NYM
        DO K=0,TNKZ
          DO I=0,NKX
            CF1(I,K,J)=(CU3(I,K,J+1)-CU3(I,K,J-1))/(2.*DYF(J))
            CF2(I,K,J)=CIKZ(K)*0.5*(CU2(I,K,J)+CU2(I,K,J+1))
          END DO
        END DO
      END DO

C Convert to Physical space
      CALL FFT_XZ_TO_PHYSICAL(CF1,F1,0,NY+1)
      CALL FFT_XZ_TO_PHYSICAL(CF2,F2,0,NY+1)

C Compute product
      DO J=2,NYM
        DO K=0,NZM
          DO I=0,NXM
            F1(I,K,J)=2.*F1(I,K,J)*F2(I,K,J)
          END DO
        END DO
      END DO

      CALL FFT_XZ_TO_FOURIER(F1,CF1,0,NY+1)

C Add to RHS term
      DO J=2,NYM
        DO K=0,TNKZ
          DO I=0,NKX
            CS1(I,K,J)=CS1(I,K,J)+CF1(I,K,J)
          END DO
        END DO
      END DO

C Now get the third of the off-diagonal terms
      DO J=2,NYM
        DO K=0,TNKZ
          DO I=0,NKX
            CF1(I,K,J)=CIKZ(K)*CU1(I,K,J)
            CF2(I,K,J)=CIKX(I)*CU3(I,K,J)
          END DO
        END DO
      END DO

      CALL FFT_XZ_TO_PHYSICAL(CF1,F1,0,NY+1)
      CALL FFT_XZ_TO_PHYSICAL(CF2,F2,0,NY+1)
      
C Compute product
      DO J=2,NYM
        DO K=0,NZM
          DO I=0,NXM
            F1(I,K,J)=2.*F1(I,K,J)*F2(I,K,J)
          END DO
        END DO
      END DO

      CALL FFT_XZ_TO_FOURIER(F1,CF1,0,NY+1)

C Add to RHS term
      DO J=2,NYM
        DO K=0,TNKZ
          DO I=0,NKX
            CS1(I,K,J)=CS1(I,K,J)+CF1(I,K,J)
          END DO
        END DO
      END DO     
     
C Finally, if the buoyancy force is active, then we need to add
C the contribution of the density to the pressure.  Note that the
C plane averaged density and the corresponding hydrostatic part of the
C pressure have been cancelled, so skip the 0,0 mode
      DO N=1,N_TH

      IF (F_TYPE .EQ. 4)THEN
       DO J=2,NYM
        DO K=0,TNKZ
          DO I=0,NKX
c            IF ((I.NE.0).or.(K.NE.0)) THEN
              CS1(I,K,J)=CS1(I,K,J)+Ratio_gr*cos(ANG_BETA)*
     &          (CTH(I,K,J+1,N)-CTH(I,K,J-1,N))/(GYF(J+1)-GYF(J-1))
c            END IF
          END DO
        END DO
       END DO

       DO J=2,NYM
        DO K=0,TNKZ
          DO I=0,NKX
            CS1(I,K,J)=CS1(I,K,J) + Ratio_gr*sin(ANG_BETA)*
     &          CIKX(I)*CTH(I,K,J,N)
          END DO
        END DO
      END DO

      ELSE

      DO J=2,NYM
        DO K=0,TNKZ
          DO I=0,NKX  
            IF ((I.NE.0).or.(K.NE.0)) THEN
              CS1(I,K,J)=CS1(I,K,J)+RI_TAU(N)*
     &          (CTH(I,K,J+1,N)-CTH(I,K,J-1,N))/(GYF(J+1)-GYF(J-1))
            END IF
          END DO
        END DO
      END DO
      ENDIF
      END DO

      
C Now, the RHS term should be stored in CS1     

C Construct the tridiagonal system in Fourier space to solve for CP
C First, zero the vectors
      DO J=0,NY+1
        DO I=0,NKX
          MATL_C(I,J)=0.d0
          MATD_C(I,J)=1.d0
          MATU_C(I,J)=0.d0
          VEC_C(I,J)=(0.,0.)
        END DO
      END DO

      DO K=0,TNKZ
        DO J=2,NYM
          DO I=0,NKX
            MATL_C(I,J)=1./(DY(J)*DYF(J))
            MATD_C(I,J)=-KX2(I)-KZ2(K)-1./(DY(J+1)*DYF(J))
     &                    -1./(DY(J)*DYF(J))
            MATU_C(I,J)=1./(DY(J+1)*DYF(J))   
            VEC_C(I,J)=-1.*CS1(I,K,J)
          END DO
        END DO
        
        IF (USE_MPI) THEN
          CALL APPLY_BC_POISSON_MPI(MATL_C,MATD_C,MATU_C,VEC_C,K)
C First, do the forward sweeps
          CALL THOMAS_FORWARD_COMPLEX_MPI(MATL_C,MATD_C,MATU_C,VEC_C
     &                                 ,NY,NX/3)
C Now, do the backwared sweeps to put the solution in VEC_C
          CALL THOMAS_BACKWARD_COMPLEX_MPI(MATL_C,MATD_C,MATU_C,VEC_C
     &                                  ,NY,NX/3)
        ELSE
C Else we are running in serial mode
C Apply BCs
        DO I=0,NKX
C Use dirichlet boundary condition at the lower wall to
C prevent the tridiagonal matrix from becomming singular for i,k=0
          IF ((I.EQ.0).AND.(K.EQ.0)) THEN
            MATD_C(I,1)=1.
            MATU_C(I,1)=0.
            VEC_C(I,1)=(0.,0.)
            MATD_C(I,NY)=-1.
            MATL_C(I,NY)=1.
            VEC_C(I,NY)=(0.,0.)
          ELSE
! Here, apply Neumann boundary conditions (dp/dz=0) at the walls
            MATD_C(I,1)=1.
            MATU_C(I,1)=-1.
            VEC_C(I,1)=(0.,0.)
            MATD_C(I,NY)=-1.
            MATL_C(I,NY)=1.
            VEC_C(I,NY)=(0.,0.)
          END IF
        END DO
C Now, solve for CP
        CALL THOMAS_COMPLEX(MATL_C,MATD_C,MATU_C,VEC_C,NY,NKX)
        END IF

        DO J=1,NY
          DO I=0,NKX
            CP(I,K,J)=VEC_C(I,J)
          END DO
        END DO
      END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_TH_CHAN
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Initialize the scalar fields
C In this subroutine, you should initialize each scalar field for the
C particular problem of interest

      INCLUDE 'header'
      INTEGER I,J,K,N
      real*8  th_b,L,t1in,s_1,s_2,A,fact_1 

      DO N=1,N_TH
       IF (CREATE_NEW_TH(N)) THEN

       IF ((IC_TYPE.eq.0).or.(IC_TYPE.eq.1).or.(IC_TYPE.eq.2)) THEN
       IF(DEV_BACK_TH) THEN
        TH(:,:,:,N)=0.0d0 ;
      
        IF(N ==1)THEN
        DO J=0,NY+1
! Linear profile with slope corresponding to upper value
          TH_BACK(J,N)=Temp_0 ;
        END DO
        ELSEIF (N == 2)THEN
         DO J=0,NY+1
! Linear profile with slope corresponding to upper value
          TH_BACK(J,N)=Sal_0 ;
         END DO
        ENDIF

       ELSE

       ENDIF 
      ELSE IF (IC_TYPE.eq.3) THEN
! Shear layer
       DO J=0,NY
         DO K=0,NZM
           DO I=0,NXM
             TH(I,K,J,N)=TANH(GYF(J)*15.d0)+1.d0
            END DO
          END DO
        END DO
      ELSE IF (IC_TYPE.eq.4) THEN
!   Inclined channel
        CURVE_FITT = .true.

        IF (CURVE_FITT) THEN
         DO J=0,NY+1
          DO K=0,NZM
           DO I=0,NXM
             TH(I,K,J,N)=0.0d0 ;
            END DO
          END DO
         END DO
        ELSE
         fact1=30.d0    ;
         fact2=0.0016d0  ;
         fact_1=2.0d0   ;
      
         L=1000.d0     ;
         A=fact1*(NU*L/(4.0*sqrt(RI_FINAL - OMEGA0*OMEGA0)))**(1.0/3.0)
         u_0= - fact2*0.75*(OMEGA0/sqrt(RI_FINAL))*
     &      (sqrt(RI_FINAL-omega0*omega0)/NU)**(1.0/3.0)*L**(2.0/3.0) ; ! sing change of u_0
 
         write(6,*)'I am here', RI_FINAL, A,U_0
         fact3=-fact_1*1.0/(2.0*A) ;
         fact4=sqrt(3.0)/(2.0*A) ;
        
         t1in =-0.0d0/omega0;

         DO J=0,NY+1
          s_1=exp(fact3*gyf(j)) ;
          s_2=sin(fact4*gyf(j)) ;
         
          DO K=0,NZM
            DO I=0,NXM
              TH(I,K,J,N)=(u_0/0.47)*s_1*s_2*t1in*sin(ANG_BETA)
     :                     *Grad_rho(N) ;
            END DO
           END DO
          END DO
         ENDIF   

        th_b = GYF(0)   
        call MPI_BCAST_REAL(th_b,1,1) 
        DO J=0,NY+1
! Linear profile with slope corresponding to upper value
            TH_BACK(J,N)=(GYF(J)-th_b)
        END DO
      END IF
      
      CALL FFT_XZ_TO_FOURIER(TH(0,0,0,n),CTH(0,0,0,n),0,NY+1)
      

      END IF
      END DO


      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_FLOW_CHAN
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INCLUDE "mpif.h"

      include 'header_mpi'

      INTEGER I,J,K
      REAL*8 RNUM1,RNUM2,RNUM3,Exp_din,damp_fact,s_1,s_2,t1
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed
      real*8 p1,p2,p3,p4,p5,p6,p7,p8,p9,p10
      character*35 FPROF
C Initialize the random number generator
      CALL RANDOM_SEED(SIZE = K)
      Allocate (seed(1:K))
      seed(1:K)=10
      CALL RANDOM_SEED(PUT = seed)

C UBULK0 and KICK should be set in input.dat

C Set the laminar velocity profile in physical space
       if (rank .eq. 0) then
       write(*,*) 'UBULK0: ',UBULK0
       endif

       IF (IC_TYPE.eq.0) then
C For closed channel flow
       DO J=0,NY
         DO K=0,NZM
           DO I=0,NXM
!             U1(I,K,J)= (3./2.)*UBULK0*(1.-GYF(J)**2.)
             U1(I,K,J)= UBULK0*(-GYF(J)**2.d0+ 2.0*LY*GYF(J))/LY**2.0
             U2(I,K,J)=0.
             U3(I,K,J)=0.
           END DO
         END DO
      END DO
      else if (IC_TYPE.eq.1) then 
C For open channel flow :
       DO K=0,NZM
         DO I=0,NXM
           DO J=1,NY
!            U1(I,K,J)=-(3./2.)*UBULK0*GYF(J)**2.+3.*UBULK0*GYF(J)
!            U1(I,K,J)=(-GYF(J)**2.d0+(NU+LY**2.d0)*GYF(J)/LY)/NU
!             if(GYF(J).le.0.0288d0) then
!             U3(I,K,J)= (-60.28d0*GYF(J)**2.d0+ 3.472d0*GYF(J))
!             else
!             U3(I,K,J)=0.05d0
!             endif
            
             U1(I,K,J)=0.0
             U2(I,K,J)=0.
             U3(I,K,J)=0.
           END DO
           U1(I,K,0)=0.
           U3(I,K,0)=0.
           U1(I,K,NY+1)=0.
           U3(I,K,NY+1)=0.
         END DO
      END DO
      else if (IC_TYPE.eq.2) then
C For Couette flow:
       DO J=0,NY
         DO K=0,NZM
           DO I=0,NXM
             U1(I,K,J)=gyf(j)
             U2(I,K,J)=0.
             U3(I,K,J)=0.
           END DO
         END DO
      END DO
      else if (IC_TYPE.eq.3) then
! Shear layer
       DO J=0,NY+1
         DO K=0,NZ+1
           DO I=0,NX+1
             U1(I,K,J)=TANH(GYF(J)*20.d0)
             U2(I,K,J)=0.d0
             U3(I,K,J)=0.d0
            END DO
          END DO
        END DO
      else if (IC_TYPE.eq.4) then
! Shear layer
       DO J=0,NY+1
         DO K=0,NZ+1
           DO I=0,NX+1
             U1(I,K,J)=0.d0
             U2(I,K,J)=0.d0
             U3(I,K,J)=0.d0
            END DO
          END DO
        END DO
      end if

C Zero the ghost cells
       IF (.NOT.USE_MPI) THEN       
       DO K=0,NZM
         DO I=0,NXM
           U1(I,K,0)=0.
           U2(I,K,0)=0.
           U3(I,K,0)=0.
           U1(I,K,NY+1)=0.
           U2(I,K,NY+1)=0.
           U3(I,K,NY+1)=0.
         END DO
      END DO
      END IF
   

      CALL FFT_XZ_TO_FOURIER(U1(0,0,-1),CU1(0,0,-1),0,NY+1)
      CALL FFT_XZ_TO_FOURIER(U2(0,0,-1),CU2(0,0,-1),0,NY+1)
      CALL FFT_XZ_TO_FOURIER(U3(0,0,-1),CU3(0,0,-1),0,NY+1)

      CALL FFT_XZ_TO_FOURIER(U1(0,0,NY+1),CU1(0,0,NY+1),0,1)
      CALL FFT_XZ_TO_FOURIER(U2(0,0,NY+1),CU2(0,0,NY+1),0,1)
      CALL FFT_XZ_TO_FOURIER(U3(0,0,NY+1),CU3(0,0,NY+1),0,1)
 
!      CALL FFT_XZ_TO_FOURIER(U1(0,0,0),CU1(0,0,0),0,NY+1)
!      CALL FFT_XZ_TO_FOURIER(U2(0,0,0),CU2(0,0,0),0,NY+1)
!      CALL FFT_XZ_TO_FOURIER(U3(0,0,0),CU3(0,0,0),0,NY+1)
     
      if (rank .eq. 0) then
      WRITE(*,*) 'KICK: ',KICK, CU1(0,0,NY)
      endif 

      DO I=1,NKX
        DO J=1,NY
          DO K=1,TNKZ
C Now, give the velocity field a random perturbation
            CALL RANDOM_NUMBER(RNUM1)
            CALL RANDOM_NUMBER(RNUM2)
            CALL RANDOM_NUMBER(RNUM3)

            IF (IC_TYPE.eq.3) THEN
C If we are initializing with a shear layer 
              CU1(I,K,J)=CU1(I,K,J)
     &             +(RNUM1-0.5)*KICK*EXP(-(GYF(J)*20.d0)**2.d0)
              CU2(I,K,J)=CU2(I,K,J)
     &             +(RNUM1-0.5)*KICK*EXP(-(GY(J)*20.d0)**2.d0)
              CU3(I,K,J)=CU3(I,K,J)
     &             +(RNUM1-0.5)*KICK*EXP(-(GYF(J)*20.d0)**2.d0)
            ELSE 
              CU1(I,K,J)=CU1(I,K,J)+(RNUM1-0.5)*KICK
              CU2(I,K,J)=CU2(I,K,J)+(RNUM2-0.5)*KICK
              CU3(I,K,J)=CU3(I,K,J)+(RNUM3-0.5)*KICK
            END IF
          END DO
          IF (TNKZ.EQ.0) THEN
! Here, In the 2d case we want to add a kick to the mean in z
            K=0         
            CALL RANDOM_NUMBER(RNUM1)
            CALL RANDOM_NUMBER(RNUM2)
            CALL RANDOM_NUMBER(RNUM3)

            IF (IC_TYPE.eq.3) THEN
              CU1(I,K,J)=CU1(I,K,J)
     &             +(RNUM1-0.5)*KICK*EXP(-(GYF(J)*20.d0)**2.d0)
              CU2(I,K,J)=CU2(I,K,J)
     &             +(RNUM1-0.5)*KICK*EXP(-(GYF(J)*20.d0)**2.d0)
              CU3(I,K,J)=CU3(I,K,J)
     &             +(RNUM1-0.5)*KICK*EXP(-(GYF(J)*20.d0)**2.d0)
            ELSE
              CU1(I,K,J)=CU1(I,K,J)+(RNUM1-0.5)*KICK
              CU2(I,K,J)=CU2(I,K,J)+(RNUM2-0.5)*KICK
              CU3(I,K,J)=CU3(I,K,J)+(RNUM3-0.5)*KICK
            END IF
          END IF 

        END DO
      END DO

      exp_din=GYF(NY+1);

      CALL MPI_BCAST_REAL_L(exp_din,1,1)
      if (rank .eq. 0)then
      write(6,*)'Exp_din',rank,exp_din 
      endif 

      IF ( (IC_TYPE .EQ. 0) .OR. (IC_TYPE .EQ. 4) ) THEN
             DO J=0,NY+1
              DO K=0,TNKZ
               DO I=0,NKX
                IF ((K .EQ. 0).AND.(I.EQ.0)) THEN
                  DAMP_FACT = 1.0
                ELSE
                  IF ((J .EQ. 0).and.(rank .eq.0)) THEN
                   DAMP_FACT = 1.0
                  ELSE
                   DAMP_FACT = exp(-20.d0*dble(GYF(J))/dble(exp_din))
                  ENDIF
                  CU1(I,K,J)=CU1(I,K,J)*DAMP_FACT
                  CU3(I,K,J)=CU3(I,K,J)*DAMP_FACT
                  CU2(I,K,J)=CU2(I,K,J)*DAMP_FACT
                ENDIF
               END DO
              END DO
             END DO

           IF (RANK .eq. NPROCS-1)THEN
           DO K =0,TNKZ
             DO I=0,NKX
              CU2(I,K,NY+1) = 0.d0
             ENDDO
            ENDDO
           ENDIF

!            DO J=1,NY
!             if(GYF(J).ge.0.01d0) then
!             CU3(0,0,J)=DCMPLX(0.03d0)!(-60.28d0*GYF(J)**2.d0+ 3.472d0*GYF(J))
!             print*,'here'
!            else
!            CU3(0,0,J)=DCMPLX(3.472*GYF(J))
!             endif

!            END DO

! TEMPORARY!!! REMOVE MEAN VELOCITY
           EX_DATA = .false.
           CURVE_FITT = .false.

          IF (CURVE_FITT) THEN

           IF (EX_DATA) THEN
 
           u_0  = -0.125 ;
           I = RANK+1
            
           FPROF ='prof/Prof_'
     &        //CHAR(MOD(I,100)/10+48)
     &        //CHAR(MOD(I,10)+48) // '.txt'

           OPEN (30,file=FPROF,form='formatted',status='old')  
           do j=-1,NY+2
            read(30,*) U3_0(J)
            U3_0(J) = u_0*U3_0(J) 
           enddo
           close(30)     
           IF (RANK .eq. 0) then
             write(6,*) 'Reading external data done !!!!!!!!!!!!' 
           endif 

           ELSE
            p1  = 0.00012588 ;
            p2  = -0.0043285 ;
            p3  = 0.062863 ;
            p4  = -0.50119 ;
            p5  = 2.3851   ;
            p6  = -6.8846  ;
            p7  = 11.635   ;
            p8  = -10.248  ;
            p9  = 3.1125   ;
            p10 = 0.64483  ;
           
            U3_0 = 0.0d0   ;

            u_0  = -0.075 ;

            do j=0,NY+1
             s_1     = gyf(j)/60.0
             U3_0(j) = p1*s_1**9 + p2*s_1**8 +
     &              p3*s_1**7 + p4*s_1**6 +
     &              p5*s_1**5 + p6*s_1**4 +
     &              p7*s_1**3 + p8*s_1**2 +
     &              p9*s_1    + p10

             U3_0(J) = U3_0(J)*(1.0-exp(-100.0*gyf(j)/65.0))*
     &              exp(-0.15*gyf(j)/65.0)
            ENDDO

            s_1 = maxval(U3_0)
            CALL MPI_SECLT_MAX(s_1)
            IF (RANK .eq. 0) then
             write(6,*) 'RANK', RANK, 'u_0', u_0, s_1
            endif

            DO j=0,NY+1
             U3_0(J) = u_0*U3_0(J)/s_1
            ENDDO

           ENDIF

           DO J=0,NY+1
            CU1(0,0,J)=U3_0(J)                 !0.d0
            CU2(0,0,J)=0.d0
            CU3(0,0,J)=0.d0
           enddo
          ELSE
           DO J=0,NY+1
            s_1=exp(fact3*gyf(j)) ;
            s_2=sin(fact4*gyf(j)) ;
            CU1(0,0,J)=(u_0/0.47)*s_1*s_2      !0.d0
            CU2(0,0,J)=0.d0
            CU3(0,0,J)=0.d0
           END DO
          ENDIF
        END IF
                
      IF (USE_MPI) THEN
        CALL GHOST_CHAN_MPI
      END IF

      CALL SAVE_STATS_CHAN(.FALSE.)

      !stop      

C Apply Boundary conditions to velocity field
      IF (USE_MPI) THEN
        CALL APPLY_BC_VEL_MPI
      ELSE
        CALL APPLY_BC_VEL_LOWER
        CALL APPLY_BC_VEL_UPPER
      END IF

C Remove the divergence of the velocity field
           DO J=1,NY
!            if(GYF(J).le.2.0d0) then
!            CU3(0,0,J)=-0.0125*GYF(J)**2+0.05*GYF(J)           
!            else
!            CU3(0,0,J)=DCMPLX(.050d0)
!             CU3(0,0,J)=-0.0125*GYF(J)**2+0.05*GYF(J)
!            endif
             CU3(0,0,J)=UBULK0*tanh(GYF(J)*69.44)

            END DO

      CALL REM_DIV_CHAN

      IF (USE_MPI) THEN
        CALL GHOST_CHAN_MPI
      END IF

C Get the pressure from the poisson equation
!      CALL POISSON_P_CHAN

      IF (USE_MPI) THEN
        CALL GHOST_CHAN_MPI
      END IF

      CALL SAVE_STATS_CHAN(.FALSE.)

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE VIS_FLOW_CHAN
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Convert to physical space and output the velocity and pressure 
C to be read for visualization

      INCLUDE 'header'
      CHARACTER*35 FNAME,FNAME_TH
      INTEGER      I, J, K, n

      FNAME='zodiac.vis'
      WRITE(6,*) 'Writing flow to ',FNAME
      OPEN(UNIT=10,FILE=FNAME,STATUS="UNKNOWN",FORM="UNFORMATTED")

      FNAME_TH='zodiac_th.vis'
      WRITE(6,*) 'Writing flow to ',FNAME_TH
      OPEN(UNIT=20,FILE=FNAME_TH,STATUS="UNKNOWN",FORM="UNFORMATTED")

C Write out grid at GYF points
      open(unit=11,file='ygrid_out.txt'
     &       ,status='unknown',form='formatted')
      write(11,111) (GYF(j),j=2,NYM)
      close(11)
111   format(64(F16.8,' '))
      open(unit=12,file='xgrid_out.txt'
     &       ,status='unknown',form='formatted')
      write(12,112) (GX(i),i=0,NXM)
      close(12)
112   format(64(F16.8,' '))
      open(unit=13,file='zgrid_out.txt'
     &      ,status='unknown',form='formatted')
      write(13,113) (GZ(k),k=0,NZM)
      close(13)
113   format(64(F16.8,' '))

      IF (NUM_PER_DIR.EQ.3) THEN
C        Enter write statment here
      ELSEIF (NUM_PER_DIR.EQ.2) THEN
C Convert to physical space
      call fft_xz_to_physical(CU1(0,0,0),U1(0,0,0),0,NY+1)
      call fft_xz_to_physical(CU2(0,0,0),U2(0,0,0),0,NY+1)
      call fft_xz_to_physical(CU3(0,0,0),U3(0,0,0),0,NY+1)
      call fft_xz_to_physical(CP,P,0,NY+1)
      do n=1,N_TH
        call fft_xz_to_physical(CTH(0,0,0,n),TH(0,0,0,n),0,NY+1)
      end do
C Interpolate the vertical velocity to GYF gridpoints
        WRITE(10) ((( 
     &    REAL(U1(I,K,J)),REAL(0.5*(U2(I,K,J)+U2(I,K,J+1)))
     &    ,REAL(U3(I,K,J)),REAL(U1(I,K,J)),
     &    REAL(0.5*(U2(I,K,J)+U2(I,K,J+1))),REAL(U3(I,K,J)),
     &    REAL(P(I,K,J))
     &    ,K=0,NZM),J=2,NYM),I=0,NXM)

        WRITE(20) ((((
     &     REAL(TH(I,K,J,n))
     &        ,n=1,N_TH),K=0,NZM),J=2,NYM), I=0,NXM)

C Output velocity field for input to LIC (Line integral convolution)
      open(61,file='lic_x.dat',form='formatted',status='unknown')
      write(61,*) NZ,NY
      open(62,file='lic_y.dat',form='formatted',status='unknown')
      write(62,*) NZ,NY

      do j=1,NY
        write(61,161) ((real(U3(30,k,j))),k=0,NZM)
        write(62,161) ((real(U2(30,k,j))),k=0,NZM)
      end do
161   format(192(F8.3))
      close(61)
      close(62)

      ELSEIF (NUM_PER_DIR.EQ.1) THEN
C        Enter write statment here
      ELSEIF (NUM_PER_DIR.EQ.0) THEN
C        Enter write statment here

      END IF
101   format(10(F16.8,' '))
      CLOSE(10)
      CLOSE(20)


C Compute the discrimenant for visualization of vortices
C Note, the velocity field will be destroyed by this calculation,
C so it is important that this is only done at the end of a simulation

C NOTE:  THIS SECTION NEEDS TO BE CHECKED, THERE MAY BE AN ERROR
C IN THE CALCULATION OF THE DISCRIMINANT

      IF ((NUM_PER_DIR.EQ.2)) THEN
C First, convert back to Fourier space
      call fft_xz_to_fourier(U1(0,0,0),CU1(0,0,0),0,NY+1)
      call fft_xz_to_fourier(U2(0,0,0),CU2(0,0,0),0,NY+1)
      call fft_xz_to_fourier(U3(0,0,0),CU3(0,0,0),0,NY+1)
C First, calculate the velocity gradient tensor at GYF points
      do j=2,NYM
        do k=0,TNKZ
          do i=0,NKX
            CA21(i,k,j)=CIKX(i)*0.5*(CU2(i,k,j+1)+CU2(i,k,j))
            CA31(i,k,j)=CIKX(i)*CU3(i,k,j)
            CA12(i,k,j)=(0.5*(CU1(i,k,j+1)+CU1(i,k,j))
     &                 - 0.5*(CU1(i,k,j)+CU1(i,k,j-1)))/DYF(j) 
            CA32(i,k,j)=(0.5*(CU3(i,k,j+1)+CU3(i,k,j))
     &                 - 0.5*(CU3(i,k,j)+CU3(i,k,j-1)))/DYF(j) 
            CA13(i,k,j)=CIKZ(k)*CU1(i,k,j)
            CA23(i,k,j)=CIKZ(k)*0.5*(CU2(i,k,j+1)+CU2(i,k,j))
C Now, the following will overwrite CUi
            CA11(i,k,j)=CIKX(i)*CU1(i,k,j)
            CA22(i,k,j)=(CU2(i,k,j+1)-CU2(i,k,j))/DYF(j)
            CA33(i,k,j)=CIKZ(k)*CU3(i,k,j)
          end do
        end do
      end do
C Transform to physicl space
      call fft_xz_to_physical(CA11,A11,0,NY+1)
      call fft_xz_to_physical(CA21,A21,0,NY+1)
      call fft_xz_to_physical(CA31,A31,0,NY+1)
      call fft_xz_to_physical(CA12,A12,0,NY+1)
      call fft_xz_to_physical(CA22,A22,0,NY+1)
      call fft_xz_to_physical(CA32,A32,0,NY+1)
      call fft_xz_to_physical(CA13,A13,0,NY+1)
      call fft_xz_to_physical(CA23,A23,0,NY+1)
      call fft_xz_to_physical(CA33,A33,0,NY+1)

C Defining S_ij=(A_ij+A_ji)/2, compute Strain_rate = S_ij S_ji
      do j=2,NYM
        do k=0,NZM
          do i=0,NXM
            Strain_rate(i,k,j)=0.5*((A12(i,k,j)+A21(i,k,j))**2.
     &                             +(A13(i,k,j)+A31(i,k,j))**2. 
     &                             +(A23(i,k,j)+A32(i,k,j))**2.)
     &                    +A11(i,k,j)**2.+A22(i,k,j)**2.+A33(i,k,j)**2. 
          end do
C Compute third invariant = -det(A)
C Overwrites A11
          do i=0,NXM
            Third_ivar(I,K,J) =
     &  - A11(I,K,J)*(A22(I,K,J)*A33(I,K,J) - A23(I,K,J)*A32(I,K,J))
     &  + A12(I,K,J)*(A21(I,K,J)*A33(I,K,J) - A23(I,K,J)*A31(I,K,J))
     &  - A13(I,K,J)*(A21(I,K,J)*A32(I,K,J) - A22(I,K,J)*A31(I,K,J))
          end do

C Defining Omega_ij=(A_ij-A_ji)/2, compute Enstrophy=-(Omega_ij Omega_ji). 
C Note that this loop overwrites A22.
          do i=0,NXM
            Enstrophy  (I,K,J) = 0.5*((A12(I,K,J) - A21(I,K,J))**2+
     &        (A13(I,K,J) - A31(I,K,J))**2+ 
     &        (A23(I,K,J) - A32(I,K,J))**2)
          end do

C Compute Second_ivar.
C Note that this loop overwrites A33.
          do i=0,NXM
            Second_ivar(I,K,J)=0.5*(Enstrophy(I,K,J)-Strain_rate(I,K,J))
          end do

C Compute Discriminant.
C Note that this loop overwrites A12.
          do i=0,NXM
            Discriminant(I,K,J) = 6.75*(Third_ivar(I,K,J))**2
     &       + (Second_ivar(I,K,J))**3
          end do
        end do 
      end do 
           
      OPEN(UNIT=20,FILE='inv.vis',STATUS='UNKNOWN',FORM='UNFORMATTED')
        WRITE(20) ((( 
     &    REAL(GX(I)),REAL(GYF(J))
     &    ,REAL(GZ(K)),REAL(Discriminant(i,k,j)),
     &    REAL(Enstrophy(i,k,j)),REAL(Strain_rate(i,k,j))
     &    ,K=0,NZM),J=2,NYM),I=0,NXM)
      CLOSE(20)

      END IF

      RETURN
      END

 


      subroutine filter_chan
C This subroutine applies a filter to the highest wavenumbers
C It should be applied to the scalars in Fourier space
C The filter used is a sharpened raised cosine filter in the horizontal
C and a fourth order implicit compact filter in the vertical, with the
C parameter alpha determining the width of the vertical filtering window
C See Lele J. comp phys?

      include 'header'

      integer I,J,K,js,je,N

! Variables for horizontal filtering
      real*8 sigma(0:NKX,0:TNKZ),sigma0

! Variables for vertical filtering
      real*8 alpha
      parameter (alpha=0.48d0)
! Parameters for a larger stencil filter
      real*8 f_a,f_b,f_c

      js=0
      je=NY+1

C Set the filtering constants for the horizontal direction
      DO i=0,NKX
       DO k=0,TNKZ
        sigma0=0.5d0*(1.d0+
     &       cos(sqrt((KX(i)*LX*1.d0/float(NX))**2.d0
     &            +(KZ(k)*LZ*1.d0/float(NZ))**2.d0)))
! Apply a sharpened raised cosine filter
        sigma(i,k)=sigma0**4.d0*(35.d0-84.d0*sigma0
     &        +70.d0*sigma0**2.d0-20.d0*sigma0**3.d0)
       END DO
      END DO

      DO N=1,N_TH
C Do the spectral filtering in the horizontal
        DO K=0,TNKZ
          DO I=0,NKX
            DO J=js+1,je-1
              CTH(I,K,J,N)=CTH(I,K,J,N)*sigma(i,k)
            END DO
          END DO
        END DO
      END DO
C Set the filtering constants
      f_a=(1.d0/8.d0)*(5.d0+6.d0*alpha)
      f_b=0.5d0*(1.d0+2.d0*alpha)
      f_c=(-1.d0/8.d0)*(1.d0-2.d0*alpha)


      DO N=1,N_TH
C First, zero the tridiagonal matrix components
      DO I=0,NKX
        DO J=0,NY+1
          MATD_C(I,J)=1.d0
          MATL_C(I,J)=0.d0
          MATU_C(I,J)=0.d0
          VEC_C(I,J)=0.d0
        END DO
      END DO


C Filter the passive scalar, TH in the vertical direction
      DO K=1,TNKZ
        DO I=1,NKX
C Construct the centered difference terms
          DO J=2,NY-1
            MATL_C(I,J)=alpha
            MATD_C(I,J)=1.d0
            MATU_C(I,J)=alpha
            VEC_C(I,J)=f_a*CTH(I,K,J,N)
     &                +(f_b/2.d0)*(CTH(I,K,J+1,N)+CTH(I,K,J-1,N))
     &                +(f_c/2.d0)*(CTH(I,K,J+2,N)+CTH(I,K,J-2,N))
          END DO
C Now, construct the equations for the boundary nodes
          J=1
            MATL_C(I,J)=0.d0
            MATD_C(I,J)=1.d0
            MATU_C(I,J)=0.d0
            VEC_C(I,J)=CTH(I,K,J,N)
          J=NY
            MATL_C(I,J)=0.d0
            MATD_C(I,J)=1.d0
            MATU_C(I,J)=0.d0
            VEC_C(I,J)=CTH(I,K,J,N)
         END DO
C Now, solve the tridiagonal system
         CALL THOMAS_COMPLEX(MATL_C,MATD_C,MATU_C,VEC_C,NY,NKX)
         DO I=1,NKX
           DO J=js+1,je-1
             CTH(I,K,J,N)=VEC_C(I,J)
           END DO
         END DO
C END DO K  
       END DO

C END DO N 
       END DO
       return
       end




      subroutine courant
! This subroutine sets the timestep based on the specified CFL number
! The subroutine should be called with the velocity in physical space

      include 'header'
      INCLUDE "mpif.h"

      include 'header_mpi'

      real*8 vel
!      real*8 dt,dt_temp
      real*8 dt_x,dt_y,dt_z
      integer i,j,k
      integer imin,jmin,kmin

! Set the initial dt to some arbitrary large number
      dt=999.d0

      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            dt_x=cfl*dx(i)/abs(U1(i,k,j))
            dt_y=cfl*dy(j)/abs(U2(i,k,j))
            dt_z=cfl*dz(k)/abs(U3(i,k,j))
            dt=min(dt,dt_x,dt_y,dt_z)
          end do
        end do
      end do
 
c      write(400,*) dt, time_step
      IF (USE_MPI) THEN
       CALL COURANT_MPI_1
      ENDIF
c      write(410,*) dt, time_step      

      if (dt.le.0) then
        write(*,*) 'Error: dt<=0 in courant'
! Set DELTA_T to some small default value
        write(*,*) dt
        DELTA_T=0.0001d0
      else if (dt.ge. delta_T_lt) then
        if (rank.eq.0) then
        write(*,*) 'WARNING: DELTA_T > ',delta_T_lt,
     &           ',value capped at ', delta_T_lt, dt
        endif
        DELTA_T=delta_T_lt
        H_BAR(1)=DELTA_T*(8.0/15.0)
        H_BAR(2)=DELTA_T*(2.0/15.0)
        H_BAR(3)=DELTA_T*(5.0/15.0) 
      else
        DELTA_T=dt
        H_BAR(1)=DELTA_T*(8.0/15.0)
        H_BAR(2)=DELTA_T*(2.0/15.0)
        H_BAR(3)=DELTA_T*(5.0/15.0)
      end if

      return
      end

