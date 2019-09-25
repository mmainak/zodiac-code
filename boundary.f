
C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE APPLY_BC_1_LOWER(MATL,MATD,MATU,VEC,K)
C----*|--.---------.---------.---------.---------.---------.---------.-|-----
      INCLUDE 'header'
      INTEGER I,K

C Bottom Wall:
      IF (U_BC_YMIN.EQ.0) THEN
C Dirichlet
        DO I=0,NXM
          MATL(I,0)=0. 
          MATD(I,0)=1.
          MATU(I,0)=0.                   
          VEC(I,0)=0.

          MATL(I,1)=0. 
          MATD(I,1)=1.
          MATU(I,1)=0.                   
          VEC(I,1)=U_BC_YMIN_C1 
        END DO
      ELSE IF (U_BC_YMIN.eq.1) THEN
C Neumann
        DO I=0,NXM
          MATL(I,0)=0.
          MATD(I,0)=-1.
          MATU(I,0)=1.
          VEC(I,0)=DY(1)*U_BC_YMIN_C1
        END DO
      ELSE IF (U_BC_YMIN.eq.3) THEN
C Wall model artificial slip boundary condition
C Neumann
        DO I=0,NXM
          MATL(I,1)=0.
          MATD(I,1)=-1.
          MATU(I,1)=1.
          VEC(I,1)=DY(2)*U_BC_LOWER_NWM(I,K)
        END DO
C This gridpoint is not used here
        DO I=0,NXM
          MATL(I,0)=0.
          MATD(I,0)=1.
          MATU(I,0)=0.
          VEC(I,0)=0.
        END DO
      ELSE
        write(*,*) 'WARNING: U_BC_LOWER is of unknown type'
      END IF

      RETURN 
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE APPLY_BC_1_LOWER_C(MATL_C,MATD_C,MATU_C,VEC_C)
C----*|--.---------.---------.---------.---------.---------.---------.-|-----
      INCLUDE 'header'
      INTEGER I

C Bottom Wall:
      IF (U_BC_YMIN.EQ.0) THEN
C Dirichlet
        DO I=0,NKX
          MATL_C(I,0)=0. 
          MATD_C(I,0)=1.
          MATU_C(I,0)=0.                   
          VEC_C(I,0)=0.

          MATL_C(I,1)=0. 
          MATD_C(I,1)=1.
          MATU_C(I,1)=0.                   
          VEC_C(I,1)=U_BC_YMIN_C1 
        END DO
      ELSE
C Neumann
        DO I=0,NKX
          MATL_C(I,0)=0.
          MATD_C(I,0)=-1.
          MATU_C(I,0)=1.
          VEC_C(I,0)=DY(1)*U_BC_YMIN_C1
        END DO
      END IF

      RETURN 
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|----
      SUBROUTINE APPLY_BC_1_UPPER(MATL,MATD,MATU,VEC,K)
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      INCLUDE 'header'
      INTEGER I,K

C Top wall
      IF (U_BC_YMAX.EQ.0) THEN
C Dirichlet
        DO I=0,NXM
          MATL(I,NY+1)=0.
          MATD(I,NY+1)=1.
          MATU(I,NY+1)=0.
          VEC(I,NY+1)=0.

          MATL(I,NY)=0.
          MATD(I,NY)=1.
          MATU(I,NY)=0.
          VEC(I,NY)=U_BC_YMAX_C1
        END DO
      ELSE IF (U_BC_YMAX.eq.1) THEN
C Neumann
        DO I=0,NXM
          MATL(I,NY)=-1.
          MATD(I,NY)=1.
          MATU(I,NY)=0.
          VEC(I,NY)=DY(NY)*U_BC_YMAX_C1
        END DO
        DO I=0,NXM
          MATL(I,NY+1)=0.
          MATD(I,NY+1)=1.
          MATU(I,NY+1)=0.
          VEC(I,NY+1)=0.
        END DO
      ELSE IF (U_BC_YMAX.EQ.3) THEN
C Use wall model artificial boundary condition
C Neumann
        DO I=0,NXM
          MATL(I,NY)=-1.
          MATD(I,NY)=1.
          MATU(I,NY)=0.
          VEC(I,NY)=DY(NY)*U_BC_UPPER_NWM(I,K)
        END DO
C This gridpoint is not used here
        DO I=0,NXM
          MATL(I,NY+1)=0.
          MATD(I,NY+1)=1.
          MATU(I,NY+1)=0.
          VEC(I,NY+1)=U1(I,K,NY)
        END DO
      END IF

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|----
      SUBROUTINE APPLY_BC_1_UPPER_C(MATL_C,MATD_C,MATU_C,VEC_C)
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      INCLUDE 'header'
      INTEGER I

C Top wall
      IF (U_BC_YMAX.EQ.0) THEN
C Dirichlet
        DO I=0,NKX
          MATL_C(I,NY+1)=0.
          MATD_C(I,NY+1)=1.
          MATU_C(I,NY+1)=0.
          VEC_C(I,NY+1)=0.

          MATL_C(I,NY)=0.
          MATD_C(I,NY)=1.
          MATU_C(I,NY)=0.
          VEC_C(I,NY)=U_BC_YMAX_C1
        END DO
      ELSE
C Neumann
        DO I=0,NKX
          MATL_C(I,NY+1)=-1.
          MATD_C(I,NY+1)=1.
          MATU_C(I,NY+1)=0.
          VEC_C(I,NY+1)=DY(NY+1)*U_BC_YMAX_C1
        END DO      
      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|---
      SUBROUTINE APPLY_BC_2_LOWER(MATL,MATD,MATU,VEC,K)
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      INCLUDE 'header'
      INTEGER I,K

C Bottom Wall:
      IF (V_BC_YMIN.EQ.0) THEN
C Dirichlet
        DO I=0,NXM
          MATL(I,1)=0.d0 
          MATD(I,1)=1.d0
          MATU(I,1)=0.d0                   
          VEC(I,1)=V_BC_YMIN_C1 

          MATL(I,2)=0.d0 
          MATD(I,2)=1.d0
          MATU(I,2)=0.d0                   
          VEC(I,2)=V_BC_YMIN_C1 
        END DO
      ELSE IF (V_BC_YMIN.EQ.1) THEN
C Neumann
        DO I=0,NXM
          MATD(I,1)=-1.d0
          MATU(I,1)=1.d0
          MATL(I,1)=0.d0
          VEC(I,1)=DYF(1)*V_BC_YMIN_C1
        END DO
      ELSE IF (V_BC_YMIN.EQ.3) THEN
C Wall model, wall located at GY(2)
        DO I=0,NXM
          MATL(I,1)=0.d0
          MATD(I,1)=1.d0
          MATU(I,1)=0.d0
          VEC(I,1)=0.d0

          MATL(I,2)=0.d0
          MATD(I,2)=1.d0
          MATU(I,2)=0.d0
          VEC(I,2)=0.d0
        END DO
      END IF

C The following is only a placeholder, this row is used for U1 and U3
      DO I=0,NXM
        MATL(I,0) = 0.
        MATD(I,0) = 1.
        MATU(I,0) = 0.
        VEC(I,0) = 0.
      END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|---
      SUBROUTINE APPLY_BC_2_LOWER_C(MATL_C,MATD_C,MATU_C,VEC_C)
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      INCLUDE 'header'
      INTEGER I

C Bottom Wall:
      IF (V_BC_YMIN.EQ.0) THEN
C Dirichlet
        DO I=0,NKX
          MATL_C(I,1)=0.d0 
          MATD_C(I,1)=1.d0
          MATU_C(I,1)=0.d0                   
          VEC_C(I,1)=V_BC_YMIN_C1 

          MATL_C(I,2)=0.d0 
          MATD_C(I,2)=1.d0
          MATU_C(I,2)=0.d0                   
          VEC_C(I,2)=V_BC_YMIN_C1 
        END DO
      ELSE IF (V_BC_YMIN.EQ.1) THEN
C Neumann
        DO I=0,NKX
          MATD_C(I,1)=-1.d0
          MATU_C(I,1)=1.d0
          MATL_C(I,1)=0.d0
          VEC_C(I,1)=DYF(1)*V_BC_YMIN_C1
        END DO
      END IF

C The following is only a placeholder, this row is used for U1 and U3
      DO I=0,NKX
        MATL_C(I,0) = 0.
        MATD_C(I,0) = 1.
        MATU_C(I,0) = 0.
        VEC_C(I,0) = 0.
      END DO

      RETURN
      END

 
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      SUBROUTINE APPLY_BC_2_UPPER(MATL,MATD,MATU,VEC,K)
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      INCLUDE 'header'
      INTEGER I,K
C Top wall
      IF (V_BC_YMAX.EQ.0) THEN
C Dirichlet
        DO I=0,NXM
          MATL(I,NY+1)=0.
          MATD(I,NY+1)=1.
          MATU(I,NY+1)=0.
          VEC(I,NY+1)=V_BC_YMAX_C1
          
          MATL(I,NY)=0.
          MATD(I,NY)=1.
          MATU(I,NY)=0.
          VEC(I,NY)=V_BC_YMAX_C1
        END DO
      ELSE IF (V_BC_YMAX.EQ.1) THEN
C Neumann
        DO I=0,NXM
          MATL(I,NY+1)=-1.
          MATD(I,NY+1)=1.
          MATU(I,NY+1)=0.
          VEC(I,NY+1)=DYF(NY)*V_BC_YMAX_C1
        END DO
      ELSE IF (V_BC_YMAX.eq.3) THEN
C Wall model, wall located at GY(NY)
        DO I=0,NXM
          MATL(I,NY+1)=0.
          MATD(I,NY+1)=1.
          MATU(I,NY+1)=0.
          VEC(I,NY+1)=0.

          MATL(I,NY)=0.
          MATD(I,NY)=1.
          MATU(I,NY)=0.
          VEC(I,NY)=0.
        END DO
      END IF      
      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|--
      SUBROUTINE APPLY_BC_2_UPPER_C(MATL_C,MATD_C,MATU_C,VEC_C)
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      INCLUDE 'header'
      INTEGER I
C Top wall
      IF (V_BC_YMAX.EQ.0) THEN
C Dirichlet
        DO I=0,NKX
          MATL_C(I,NY+1)=0.
          MATD_C(I,NY+1)=1.
          MATU_C(I,NY+1)=0.
          VEC_C(I,NY+1)=V_BC_YMAX_C1
          
          MATL_C(I,NY)=0.
          MATD_C(I,NY)=1.
          MATU_C(I,NY)=0.
          VEC_C(I,NY)=V_BC_YMAX_C1
        END DO
      ELSE IF (V_BC_YMAX.EQ.1) THEN
C Neumann
        DO I=0,NKX
          MATL_C(I,NY+1)=-1.
          MATD_C(I,NY+1)=1.
          VEC_C(I,NY+1)=DYF(NY)*V_BC_YMAX_C1
        END DO      
      END IF
      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|--
      SUBROUTINE APPLY_BC_3_LOWER(MATL,MATD,MATU,VEC,K)
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      INCLUDE 'header'
      INTEGER I,K

C Bottom Wall:
      IF (W_BC_YMIN.EQ.0) THEN
C Dirichlet
        DO I=0,NXM
          MATL(I,0)=0. 
          MATD(I,0)=1.
          MATU(I,0)=0.                   
          VEC(I,0)=0.

          MATL(I,1)=0. 
          MATD(I,1)=1.
          MATU(I,1)=0.                   
          VEC(I,1)=W_BC_YMIN_C1
        END DO
      ELSE IF (W_BC_YMIN.EQ.1) THEN
C Neumann
        DO I=0,NXM
          MATL(I,0)=0.
          MATD(I,0)=-1.
          MATU(I,0)=1.
          VEC(I,0)=DY(1)*W_BC_YMIN_C1
        END DO
      ELSE IF (W_BC_YMIN.EQ.3) THEN
C Wall model, artificial slip boundary condition
C Apply at GY(2), which is now considered the wall location
C Neumann
        DO I=0,NXM
          MATL(I,1)=0.
          MATD(I,1)=-1.
          MATU(I,1)=1.
          VEC(I,1)=DY(2)*W_BC_LOWER_NWM(I,K)
        END DO
! This gridpoint is not used for the near wall model
        DO I=0,NXM
          MATL(I,0)=0.
          MATD(I,0)=1.
          MATU(I,0)=0.
          VEC(I,0)=0.
        END DO
      END IF

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|--
      SUBROUTINE APPLY_BC_3_LOWER_C(MATL_C,MATD_C,MATU_C,VEC_C)
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      INCLUDE 'header'
      INTEGER I

C Bottom Wall:
      IF (W_BC_YMIN.EQ.0) THEN
C Dirichlet
        DO I=0,NKX
          MATL_C(I,0)=0. 
          MATD_C(I,0)=1.
          MATU_C(I,0)=0.                   
          VEC_C(I,0)=0.

          MATL_C(I,1)=0. 
          MATD_C(I,1)=1.
          MATU_C(I,1)=0.                   
          VEC_C(I,1)=W_BC_YMIN_C1
        END DO
      ELSE
C Neumann
        DO I=0,NKX
          MATL_C(I,0)=0.
          MATD_C(I,0)=-1.
          MATU_C(I,0)=1.
          VEC_C(I,0)=DY(1)*W_BC_YMIN_C1
        END DO
      END IF

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|--
      SUBROUTINE APPLY_BC_3_UPPER(MATL,MATD,MATU,VEC,K)
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      INCLUDE 'header'
      INTEGER I,K

C Top wall
      IF (W_BC_YMAX.EQ.0) THEN
C Dirichlet
        DO I=0,NXM
          MATL(I,NY+1)=0.
          MATD(I,NY+1)=1.
          MATU(I,NY+1)=0.
          VEC(I,NY+1)=0.

          MATL(I,NY)=0.
          MATD(I,NY)=1.
          MATU(I,NY)=0.
          VEC(I,NY)=W_BC_YMAX_C1
        END DO
      ELSE IF (W_BC_YMAX.EQ.1) THEN
C Neumann
        DO I=0,NXM
          MATL(I,NY)=-1.
          MATD(I,NY)=1.
          MATU(I,NY)=0.
          VEC(I,NY)=DY(NY)*W_BC_YMAX_C1
        END DO
        DO I=0,NXM
          MATL(I,NY+1)=0.
          MATD(I,NY+1)=1.
          MATU(I,NY+1)=0.
          VEC(I,NY+1)=0.
        END DO
      ELSE IF (W_BC_YMAX.EQ.3) THEN
C Wall model artificial slip boundary condition
C Neumann
        DO I=0,NXM
          MATL(I,NY)=-1.
          MATD(I,NY)=1.
          MATU(I,NY)=0.
          VEC(I,NY)=DY(NY)*W_BC_UPPER_NWM(I,K)
        END DO
C This gridpoint is not used here
        DO I=0,NXM
          MATL(I,NY+1)=0.
          MATD(I,NY+1)=1.
          MATU(I,NY+1)=0.
          VEC(I,NY+1)=0.
        END DO

      END IF

      RETURN 
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|--
      SUBROUTINE APPLY_BC_3_UPPER_C(MATL_C,MATD_C,MATU_C,VEC_C)
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      INCLUDE 'header'
      INTEGER I

C Top wall
      IF (W_BC_YMAX.EQ.0) THEN
C Dirichlet
        DO I=0,NKX
          MATL_C(I,NY+1)=0.
          MATD_C(I,NY+1)=1.
          MATU_C(I,NY+1)=0.
          VEC_C(I,NY+1)=0.

          MATL_C(I,NY)=0.
          MATD_C(I,NY)=1.
          MATU_C(I,NY)=0.
          VEC_C(I,NY)=W_BC_YMAX_C1
        END DO
      ELSE
C Neumann
        DO I=0,NKX
          MATL_C(I,NY+1)=-1.
          MATD_C(I,NY+1)=1.
          MATU_C(I,NY+1)=0.
          VEC_C(I,NY+1)=DY(NY+1)*W_BC_YMAX_C1
        END DO      
      END IF

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|--
      subroutine APPLY_BC_TH_LOWER(MATL,MATD,MATU,VEC,K,N)
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      include 'header'
      integer i,k, N
! Bottom Wall:
      if (TH_BC_YMIN(N).eq.0) then
! Dirichlet
        do i=0,NXM
          MATL(i,0)=0. 
          MATD(i,0)=1.
          MATU(i,0)=0.                   
          VEC(i,0)=TH_BC_YMIN_C1(N)

          MATL(i,1)=0. 
          MATD(i,1)=1.
          MATU(i,1)=0.                   
          VEC(i,1)=TH_BC_YMIN_C1(N)
        end do
      elseif (TH_BC_YMIN(N).eq.3) then
! Melting boundary condition 
         IF (N==1) THEN
          do i=0,NXM
          MATL(i,0)=0.
          MATD(i,0)=1.
          MATU(i,0)=0.
          VEC(i,0) =TH_BC_YMIN_C1(N)*((TH(I,K,0,2) + TH_BACK(0,2))* 
     &                   m_FP - TH_BACK(0,1))

          MATL(i,1)=0.
          MATD(i,1)=1.
          MATU(i,1)=0.
          VEC(i,1) =TH_BC_YMIN_C1(N)*((TH(I,K,1,2) + TH_BACK(1,2))*
     &                   m_FP - TH_BACK(1,1))
          end do
         ELSE
          do i=0,NXM
          MATL(i,0)=0.
          MATD(i,0)=-1.
          MATU(i,0)=1.
          VEC(i,0)=DY(1)*TH_BC_YMIN_C1(N)*
     &               (5.0d0*PR(2)/PR(1))*(C_sp_heat/L_heat)* 
     &               (TH(I,K,1,1) + TH_BACK(1,1))*dthdz_mean(1)/m_FP
          enddo

         ENDIF   
      else
! Neumann
! NOTE: BC enforced at GY(2)
       IF (F_TYPE .EQ. 4) THEN
         do i=0,NXM
          MATL(i,1)=0.
          MATD(i,1)=-1.
          MATU(i,1)=1.
          VEC(i,1)=DY(2)*TH_BC_YMIN_C1(N)*cos(ANG_BETA)*Grad_rho(N)
        end do

        do i=0,NXM
          MATL(i,0)=0.
          MATD(i,0)=-1.
          MATU(i,0)=1.
          VEC(i,0)=DY(1)*TH_BC_YMIN_C1(N)*cos(ANG_BETA)*GRAD_rho(N)
        end do

       ELSE 
        do i=0,NXM
          MATL(i,1)=0.
          MATD(i,1)=-1.
          MATU(i,1)=1.
          VEC(i,1)=DY(2)*TH_BC_YMIN_C1(N)
        end do

        do i=0,NXM
          MATL(i,0)=0.
          MATD(i,0)=-1.
          MATU(i,0)=1.
          VEC(i,0)=DY(1)*TH_BC_YMIN_C1(N)
        end do
       ENDIF 

      end if
      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|--
      subroutine APPLY_BC_TH_LOWER_C(MATL_C,MATD_C,MATU_C,VEC_C,N)
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      include 'header'
      integer i,N
! Bottom Wall:
      if (TH_BC_YMIN(N).eq.0) then
! Dirichlet
        do i=0,NKX
          MATL_C(i,0)=0. 
          MATD_C(i,0)=1.
          MATU_C(i,0)=0.                   
          VEC_C(i,0)=0.

          MATL_C(i,1)=0. 
          MATD_C(i,1)=1.
          MATU_C(i,1)=0.                   
          VEC_C(i,1)=TH_BC_YMIN_C1(N)
        end do
      else
! Neumann
! NOTE: BC enforced at GY(1)
        do i=0,NXM
          MATL(i,0)=0.
          MATD(i,0)=-1.
          MATU(i,0)=1.
          VEC(i,0)=DY(1)*TH_BC_YMIN_C1(N)
        end do


      end if
      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|--
      subroutine APPLY_BC_TH_UPPER(MATL,MATD,MATU,VEC,N)
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      include 'header'
      integer i,N
! Top wall
      if (TH_BC_YMAX(N).eq.0) then
! Dirichlet
        do i=0,NXM
          MATL(i,NY+1)=0.
          MATD(i,NY+1)=1.
          MATU(i,NY+1)=0.
          VEC(i,NY+1)=TH_BC_YMAX_C1(N)

          MATL(i,NY)=0.
          MATD(i,NY)=1.
          MATU(i,NY)=0.
          VEC(i,NY)=TH_BC_YMAX_C1(N)
        end do
      else
! Neumann
! NOTE: BC enforced at GY(NY+1)
        do i=0,NXM
          MATL(i,NY)=-1.
          MATD(i,NY)=1.
          MATU(i,NY)=0.
          VEC(i,NY)=DY(NY)*TH_BC_YMAX_C1(N)
        end do
        do i=0,NXM
          MATL(i,NY+1)=-1.
          MATD(i,NY+1)=1.
          MATU(i,NY+1)=0.
          VEC(i,NY+1)=DY(NY+1)*TH_BC_YMAX_C1(N)
        end do


      end if
      return
      end

C----*|--.---------.---------.---------.---------.---------.---------.-|--
      subroutine APPLY_BC_TH_UPPER_C(MATL_C,MATD_C,MATU_C,VEC_C,N)
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      include 'header'
      integer i,N
! Top wall
      if (TH_BC_YMAX(N).eq.0) then
! Dirichlet
        do i=0,NKX
          MATL_C(i,NY+1)=0.
          MATD_C(i,NY+1)=1.
          MATU_C(i,NY+1)=0.
          VEC_C(i,NY+1)=0.

          MATL_C(i,NY)=0.
          MATD_C(i,NY)=1.
          MATU_C(i,NY)=0.
          VEC_C(i,NY)=TH_BC_YMAX_C1(N)
        end do
      else
! Neumann
! NOTE: BC enforced at GY(NY+1)
        do i=0,NKX
          MATL_C(i,NY+1)=-1.
          MATD_C(i,NY+1)=1.
          MATU_C(i,NY+1)=0.
          VEC_C(i,NY+1)=DY(NY+1)*TH_BC_YMAX_C1(N)
        end do      
      end if
      return
      end

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE THOMAS_REAL(A,B,C,G,NY,NX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Uses the Thomas algorithm to solve Ax=b for tridiagonal A
C The RHS vector and solution are real
C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
C Returns solution in x
C The indexing should be done by ROW, ie.
C [ b1  c1   0   0   0 ...
C [ a2  b2  c2   0   0 ...
C [  0  a3  b3   c3  0 ...

      INTEGER I, J, NX, NY
      REAL*8 A(0:NX,0:NY), B(0:NX,0:NY), C(0:NX,0:NY), G(0:NX,0:NY)

      DO J=0,NY-1
        DO I=0,NX
          A(I,J+1)=-A(I,J+1)/B(I,J)
          B(I,J+1)=B(I,J+1)+A(I,J+1)*C(I,J)
          G(I,J+1)=G(I,J+1)+A(I,J+1)*G(I,J)
        END DO
      END DO
      DO I=0,NX
        G(I,NY)=G(I,NY)/B(I,NY)
      END DO
      DO J=NY-1,0,-1
        DO I=0,NX
          G(I,J)=(G(I,J)-C(I,J)*G(I,J+1))/B(I,J)
        END DO
      END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|    
      SUBROUTINE THOMAS_COMPLEX(A,B,C,G,NY,NX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

C Uses the Thomas algorithm to solve Ax=b for tridiagonal A
C The RHS vector and solution is complex
C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
C Returns solution in x
C The indexing should be done by ROW, ie.
C [ b1  c1   0   0   0 ...
C [ a2  b2  c2   0   0 ...
C [  0  a3  b3   c3  0 ...

      INTEGER I, J, NY, NX
      REAL*8 A(0:NX,0:NY), B(0:NX,0:NY), C(0:NX,0:NY)
      COMPLEX*16 G(0:NX,0:NY)

      DO J=0,NY-1
        DO I=0,NX
          A(I,J+1)=-A(I,J+1)/B(I,J)
          B(I,J+1)=B(I,J+1)+A(I,J+1)*C(I,J)
          G(I,J+1)=G(I,J+1)+A(I,J+1)*G(I,J)
        END DO
      END DO
      DO I=0,NX
        G(I,NY)=G(I,NY)/B(I,NY)
      END DO
      DO I=0,NX
        DO J=NY-1,0,-1
          G(I,J)=(G(I,J)-C(I,J)*G(I,J+1))/B(I,J)
        END DO
      END DO

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|--
      SUBROUTINE APPLY_BC_VEL_LOWER
C----*|--.---------.---------.---------.---------.---------.---------.-|--
C This subroutine is called after initializing the flow
C It sets the appropriate boundary conditions including ghost cell values
C  on the velocity field in Fourier space
      INCLUDE 'header'
      INTEGER I,K      

C Now, apply the boundary conditions depending on the type specified 
      IF (U_BC_YMIN.EQ.0) THEN
C Dirichlet 
C Start with zero
         DO K=0,TNKZ
           DO I=0,NKX
             CU1(I,K,1)=0.d0
           END DO
         END DO
C Now, set only the mean
         CU1(0,0,1)=U_BC_YMIN_C1
C Ghost cell not used
         CU1(0,0,0)=0.d0
         CU1(0,0,-1)=0.d0
      ELSE IF (U_BC_YMIN.EQ.1) THEN
C Neumann
         DO K=0,TNKZ
           DO I=0,NKX
             CU1(I,K,0)=CU1(I,K,1)-DY(1)*U_BC_YMIN_C1 
           END DO
         END DO
      ELSE
         STOP 'Error: U_BC_YMIN must be 0, or 1'
      END IF

      IF (W_BC_YMIN.EQ.0) THEN
C Dirichlet
C Start with zero
         DO K=0,TNKZ
           DO I=0,NKX
             CU3(I,K,1)=0.d0
           END DO
         END DO
C Now, set only the mean
         CU3(0,0,1)=W_BC_YMIN_C1
C Ghost cell not used
         CU3(0,0,0)=0.d0
         CU3(0,0,-1)=0.d0
      ELSE IF (W_BC_YMIN.EQ.1) THEN
C Neumann
         DO K=0,TNKZ
           DO I=0,NKX
             CU3(I,K,0)=CU3(I,K,1)-DY(1)*W_BC_YMIN_C1 
           END DO
         END DO
      ELSE
         STOP 'Error: W_BC_YMIN must be 0, or 1' 
      END IF

      IF (V_BC_YMIN.EQ.0) THEN
C Dirichlet
C Set the vertical velocity at GYF(1) (halfway between GY(2) and GY(1))
         DO K=0,TNKZ
           DO I=0,NKX
             CU2(I,K,1)=2.d0*V_BC_YMIN_C1-CU2(I,K,2)
             CU2(I,K,0)=2.d0*CU2(I,K,1)-CU2(I,K,2) 
           END DO
         END DO
      ELSE IF (V_BC_YMIN.EQ.1) THEN
C Neumann
         DO K=0,TNKZ
           DO I=0,NKX
             CU2(I,K,1)=CU2(I,K,2)-DYF(1)*V_BC_YMIN_C1 
           END DO
         END DO
      ELSE IF (V_BC_YMIN.EQ.2) THEN
C Upstream-travelling wave proposed by Speyer/Kim
C (initialize as zero)
        CU2(0,0,1)=-CU2(0,0,2)
      ELSE
         STOP 'Error: V_BC_YMIN must be 0, 1, or 2'
      END IF

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|--
      SUBROUTINE APPLY_BC_VEL_UPPER
C----*|--.---------.---------.---------.---------.---------.---------.-|--
C This subroutine is called after initializing the flow
C It sets the appropriate boundary conditions including ghost cell values
C  on the velocity field in Fourier space
      INCLUDE 'header'
      INTEGER I,K      

! Now, apply boundary conditions to the top of the domain
      IF (U_BC_YMAX.EQ.0) THEN
C Dirichlet 
C Start with zero
         DO K=0,TNKZ
           DO I=0,NKX
             CU1(I,K,NY)=0.d0
           END DO
         END DO
C Now, set only the mean
         CU1(0,0,NY)=U_BC_YMAX_C1
C Ghost cell not used
         CU1(0,0,NY+1)=0.d0
         CU1(0,0,NY+2)=0.d0
      ELSE IF (U_BC_YMAX.EQ.1) THEN
C Neumann
         DO K=0,TNKZ
           DO I=0,NKX
             CU1(I,K,NY+1)=CU1(I,K,NY)+DY(NY+1)*U_BC_YMAX_C1
           END DO
         END DO
      ELSE
         STOP 'Error: U_BC_YMAX must be 0, or 1'
      END IF

      IF (W_BC_YMAX.EQ.0) THEN
C Dirichlet
C Start with zero
         DO K=0,TNKZ
           DO I=0,NKX
             CU3(I,K,NY)=0.d0
           END DO
         END DO
C Now, set only the mean
         CU3(0,0,NY)=W_BC_YMAX_C1
C Ghost cell not used
         CU3(0,0,NY+1)=0.d0
         CU3(0,0,NY+2)=0.d0
      ELSE IF (W_BC_YMAX.EQ.1) THEN
C Neumann
         DO K=0,TNKZ
           DO I=0,NKX
             CU3(I,K,NY+1)=CU3(I,K,NY)+DY(NY+1)*W_BC_YMAX_C1
           END DO
         END DO
      ELSE
        STOP 'Error: W_BC_YMAX must be 0, or 1'
      END IF

      IF (V_BC_YMAX.EQ.0) THEN
C Dirichlet
C Set the vertical velocity at GYF(NY) (halfway between GY(NY) and GY(NY+1))
         DO K=0,TNKZ
           DO I=0,NKX
             CU2(I,K,NY+1)=2.d0*V_BC_YMAX_C1-CU2(I,K,NY)
             CU2(I,K,NY+2)=2.d0*CU2(I,K,NY+1)-CU2(I,K,NY)
           END DO
         END DO
      ELSE IF (V_BC_YMAX.EQ.1) THEN
C Neumann
         DO K=0,TNKZ
           DO I=0,NKX
             CU2(I,K,NY+1)=CU2(I,K,NY)+DYF(NY)*V_BC_YMAX_C1
           END DO
         END DO
      ELSE IF (V_BC_YMAX.EQ.2) THEN
C Upstream-travelling wave proposed by Speyer/Kim
C (initialize as zero gradient)
        CU2(0,0,NY+1)=-CU2(0,0,NY)
      ELSE
         STOP 'Error: V_BC_YMAX must be 0, 1, or 2'
      END IF

      RETURN
      END
