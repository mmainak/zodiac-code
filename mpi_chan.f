      SUBROUTINE GHOST_CHAN_MPI_working
! This subroutine is part of the MPI package for the channel flow
! Diablo package.
! Here, we define a set of ghost cells on each process
! the ghost cells contain information from the neighboring nodes
! and allow us to compute finite differences over the local gridpoints.
! We need to update the contents of the ghost cells at the start of
! each Runge-Kutta substep

      include 'header'
      INCLUDE "mpif.h"
      include 'header_mpi'

      integer i,j,k,N

! Define the arrays that will be used for data packing.  This makes the
! communication between processes more efficient by only requiring one
! send and recieve.
! The communication will be done in Fourier space, so these arrays should
! be complex arrays to match the velocity
! The size of the buffer arrys is 0:NKX,0:TNKZ,# of variables
      COMPLEX*16 OCPACK(0:NX/3,0:2*NZ/3,4+N_TH)
      COMPLEX*16 ICPACK(0:NX/3,0:2*NZ/3,4+N_TH)

! If we are using more than one processor, then we need to pass data
      IF (NPROCS.gt.1) THEN

! First, Pass data up the chain to higher ranked processes

      IF (RANK.eq.0) THEN
! If we are the lowest ranked process, then we don't need to recieve
! data at the lower ghost cells, these will be filled with boundary
! condition information
        DO K=0,TNKZ
          DO I=0,NKX
            OCPACK(I,K,1)=CU1(I,K,NY)
            OCPACK(I,K,2)=CU2(I,K,NY)
            OCPACK(I,K,3)=CU3(I,K,NY)
            OCPACK(I,K,4)=CP(I,K,NY)
            DO N=1,N_TH
              OCPACK(I,K,4+N)=CTH(I,K,NY,N)
            END DO
          END DO
        END DO
! Now, we have packed the data into a compact array, pass the data up
        CALL MPI_SEND(OCPACK,(4+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK+1,1,MPI_COMM_WORLD,IERROR)

! End if RANK=0
      ELSE IF (RANK.lt.NPROCS-1) THEN
! Here, we are one of the middle processes and we need to pass data
! up and recieve data from below
        DO K=0,TNKZ
          DO I=0,NKX
            OCPACK(I,K,1)=CU1(I,K,NY)
            OCPACK(I,K,2)=CU2(I,K,NY)
            OCPACK(I,K,3)=CU3(I,K,NY)
            OCPACK(I,K,4)=CP(I,K,NY)
            DO N=1,N_TH
              OCPACK(I,K,4+N)=CTH(I,K,NY,N)
            END DO
          END DO
        END DO
! Use MPI_SENDRECV since we need to recieve and send data
        CALL MPI_SEND(OCPACK,(4+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK+1,1,MPI_COMM_WORLD,IERROR)

        CALL MPI_RECV(ICPACK,(4+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK-1,1,MPI_COMM_WORLD,STATUS,IERROR)
! Now, unpack the data that we have recieved
        DO K=0,TNKZ
          DO I=0,NKX
            CU1(I,K,1)=ICPACK(I,K,1)
            CU2(I,K,1)=ICPACK(I,K,2)
            CU3(I,K,1)=ICPACK(I,K,3)
            CP(I,K,1)=ICPACK(I,K,4)
            DO N=1,N_TH
              CTH(I,K,1,N)=ICPACK(I,K,4+N)
            END DO
          END DO
        END DO

      ELSE
! Otherwise, we must be the uppermost process with RANK=NPROCS-1
! Here, we need to recieve data from below, but don't need to send data up
        CALL MPI_RECV(ICPACK,(4+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK-1,1,MPI_COMM_WORLD,STATUS,IERROR)
! Unpack the data that we have recieved
        DO K=0,TNKZ
          DO I=0,NKX
            CU1(I,K,1)=ICPACK(I,K,1)
            CU2(I,K,1)=ICPACK(I,K,2)
            CU3(I,K,1)=ICPACK(I,K,3)
            CP(I,K,1)=ICPACK(I,K,4)
            DO N=1,N_TH
              CTH(I,K,1,N)=ICPACK(I,K,4+N)
            END DO
          END DO
        END DO
      END IF

! AT this point we have passed data up the chain

! Now, we need to pass data down the chain of processes
      IF (RANK.eq.NPROCS-1) THEN
! If we are the higest ranked process, then we don't need to recieve
! data at the upper ghost cells, these will be filled with boundary
! condition information
        DO K=0,TNKZ
          DO I=0,NKX
            OCPACK(I,K,1)=CU1(I,K,2)
            OCPACK(I,K,2)=CU2(I,K,2)
            OCPACK(I,K,3)=CU3(I,K,2)
            OCPACK(I,K,4)=CP(I,K,2)
            DO N=1,N_TH
              OCPACK(I,K,4+N)=CTH(I,K,2,N)
            END DO
          END DO
        END DO
! Now, we have packed the data into a compact array, pass the data up
        CALL MPI_SEND(OCPACK,(4+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK-1,3,MPI_COMM_WORLD,IERROR)
      ELSE IF (RANK.GT.0) THEN
! Here, we are one of the middle processes and we need to pass data
! down and recieve data from above us
        DO K=0,TNKZ
          DO I=0,NKX
            OCPACK(I,K,1)=CU1(I,K,2)
            OCPACK(I,K,2)=CU2(I,K,2)
            OCPACK(I,K,3)=CU3(I,K,2)
            OCPACK(I,K,4)=CP(I,K,2)
            DO N=1,N_TH
              OCPACK(I,K,4+N)=CTH(I,K,2,N)
            END DO
          END DO
        END DO

        CALL MPI_SEND(OCPACK,(4+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK-1,3,MPI_COMM_WORLD,IERROR)

        CALL MPI_RECV(ICPACK,(4+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK+1,3,MPI_COMM_WORLD,STATUS,IERROR)
! Now, unpack the data that we have recieved
        DO K=0,TNKZ
          DO I=0,NKX
            CU1(I,K,NY+1)=ICPACK(I,K,1)
            CU2(I,K,NY+1)=ICPACK(I,K,2)
            CU3(I,K,NY+1)=ICPACK(I,K,3)
            CP(I,K,NY+1)=ICPACK(I,K,4)
            DO N=1,N_TH
              CTH(I,K,NY+1,N)=ICPACK(I,K,4+N)
            END DO
          END DO
        END DO
      ELSE
! Here, we must be the lowest process (RANK=0) and we need to recieve
! data from above
        CALL MPI_RECV(ICPACK,(4+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK+1,3,MPI_COMM_WORLD,STATUS,IERROR)
! Unpack the data that we have recieved
        DO K=0,TNKZ
          DO I=0,NKX
            CU1(I,K,NY+1)=ICPACK(I,K,1)
            CU2(I,K,NY+1)=ICPACK(I,K,2)
            CU3(I,K,NY+1)=ICPACK(I,K,3)
            CP(I,K,NY+1)=ICPACK(I,K,4)
            DO N=1,N_TH
              CTH(I,K,NY+1,N)=ICPACK(I,K,4+N)
            END DO
          END DO
        END DO
      END IF

      END IF

      RETURN
      END


      SUBROUTINE GHOST_CHAN_MPI
! This subroutine is part of the MPI package for the channel flow
! Diablo package.
! Here, we define a set of ghost cells on each process
! the ghost cells contain information from the neighboring nodes
! and allow us to compute finite differences over the local gridpoints.
! We need to update the contents of the ghost cells at the start of
! each Runge-Kutta substep

      include 'header'
      INCLUDE "mpif.h"
      include 'header_mpi'

      integer i,j,k,N

! Define the arrays that will be used for data packing.  This makes the
! communication between processes more efficient by only requiring one
! send and recieve.
! The communication will be done in Fourier space, so these arrays should
! be complex arrays to match the velocity
! The size of the buffer arrys is 0:NKX,0:TNKZ,# of variables
      COMPLEX*16 OCPACK(0:NX/3,0:2*NZ/3,5+N_TH)
      COMPLEX*16 ICPACK(0:NX/3,0:2*NZ/3,5+N_TH)

! If we are using more than one processor, then we need to pass data
      IF (NPROCS.gt.1) THEN

! First, Pass data up the chain to higher ranked processes

      IF (RANK.eq.0) THEN
! If we are the lowest ranked process, then we don't need to recieve
! data at the lower ghost cells, these will be filled with boundary
! condition information
        DO K=0,TNKZ
          DO I=0,NKX
            OCPACK(I,K,1)=CU1(I,K,NY)
            OCPACK(I,K,2)=CU2(I,K,NY)
            OCPACK(I,K,3)=CU3(I,K,NY)
            OCPACK(I,K,4)=CP(I,K,NY)
            DO N=1,N_TH
              OCPACK(I,K,4+N)=CTH(I,K,NY,N)
            END DO
c            OCPACK(i,k,5+N_TH)=cmplx(gx(i),gz(k))
c            write(410,*)gx(i),gz(k),ocpack(i,k,5+N_TH) 
          END DO
        END DO
        OCPACK(0,0,5+N_TH)=0.
! Now, we have packed the data into a compact array, pass the data up
        CALL MPI_SEND(OCPACK,(5+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK+1,1,MPI_COMM_WORLD,IERROR)

! End if RANK=0
      ELSE IF (RANK.lt.NPROCS-1) THEN
! Here, we are one of the middle processes and we need to pass data
! up and recieve data from below
        DO K=0,TNKZ
          DO I=0,NKX
            OCPACK(I,K,1)=CU1(I,K,NY)
            OCPACK(I,K,2)=CU2(I,K,NY)
            OCPACK(I,K,3)=CU3(I,K,NY)
            OCPACK(I,K,4)=CP(I,K,NY)
            DO N=1,N_TH
              OCPACK(I,K,4+N)=CTH(I,K,NY,N)
            END DO
          END DO
        END DO
        OCPACK(0,0,5+N_TH)=0.
! Use MPI_SENDRECV since we need to recieve and send data
        CALL MPI_SEND(OCPACK,(5+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK+1,1,MPI_COMM_WORLD,IERROR)

        CALL MPI_RECV(ICPACK,(5+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK-1,1,MPI_COMM_WORLD,STATUS,IERROR)
! Now, unpack the data that we have recieved
        DO K=0,TNKZ
          DO I=0,NKX
            CU1(I,K,1)=ICPACK(I,K,1)
            CU2(I,K,1)=ICPACK(I,K,2)
            CU3(I,K,1)=ICPACK(I,K,3)
            CP(I,K,1)=ICPACK(I,K,4)
            DO N=1,N_TH
              CTH(I,K,1,N)=ICPACK(I,K,4+N)
            END DO
c            write(410,*)gx(i),gz(k),icpack(i,k,5+N_TH) 
          END DO
        END DO

      ELSE
! Otherwise, we must be the uppermost process with RANK=NPROCS-1
! Here, we need to recieve data from below, but don't need to send data up
        CALL MPI_RECV(ICPACK,(5+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK-1,1,MPI_COMM_WORLD,STATUS,IERROR)
! Unpack the data that we have recieved
        DO K=0,TNKZ
          DO I=0,NKX
            CU1(I,K,1)=ICPACK(I,K,1)
            CU2(I,K,1)=ICPACK(I,K,2)
            CU3(I,K,1)=ICPACK(I,K,3)
            CP(I,K,1)=ICPACK(I,K,4)
            DO N=1,N_TH
              CTH(I,K,1,N)=ICPACK(I,K,4+N)
            END DO
          END DO
        END DO
      END IF

! AT this point we have passed data up the chain

! Now, we need to pass data down the chain of processes
      IF (RANK.eq.NPROCS-1) THEN
! If we are the higest ranked process, then we don't need to recieve
! data at the upper ghost cells, these will be filled with boundary
! condition information
        DO K=0,TNKZ
          DO I=0,NKX
            OCPACK(I,K,1)=CU1(I,K,2)
            OCPACK(I,K,2)=CU2(I,K,2)
            OCPACK(I,K,3)=CU3(I,K,2)
            OCPACK(I,K,4)=CP(I,K,2)
            DO N=1,N_TH
              OCPACK(I,K,4+N)=CTH(I,K,2,N)
            END DO
          END DO
        END DO
        OCPACK(0,0,5+N_TH)=0.
! Now, we have packed the data into a compact array, pass the data up
        CALL MPI_SEND(OCPACK,(5+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK-1,3,MPI_COMM_WORLD,IERROR)
      ELSE IF (RANK.GT.0) THEN
! Here, we are one of the middle processes and we need to pass data
! down and recieve data from above us
        DO K=0,TNKZ
          DO I=0,NKX
            OCPACK(I,K,1)=CU1(I,K,2)
            OCPACK(I,K,2)=CU2(I,K,2)
            OCPACK(I,K,3)=CU3(I,K,2)
            OCPACK(I,K,4)=CP(I,K,2)
            DO N=1,N_TH
              OCPACK(I,K,4+N)=CTH(I,K,2,N)
            END DO
          END DO
        END DO
        OCPACK(0,0,5+N_TH)=0.
        CALL MPI_SEND(OCPACK,(5+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK-1,3,MPI_COMM_WORLD,IERROR)

        CALL MPI_RECV(ICPACK,(5+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK+1,3,MPI_COMM_WORLD,STATUS,IERROR)
! Now, unpack the data that we have recieved
        DO K=0,TNKZ
          DO I=0,NKX
            CU1(I,K,NY+1)=ICPACK(I,K,1)
            CU2(I,K,NY+1)=ICPACK(I,K,2)
            CU3(I,K,NY+1)=ICPACK(I,K,3)
            CP(I,K,NY+1)=ICPACK(I,K,4)
            DO N=1,N_TH
              CTH(I,K,NY+1,N)=ICPACK(I,K,4+N)
            END DO
          END DO
        END DO
      ELSE
! Here, we must be the lowest process (RANK=0) and we need to recieve
! data from above
        CALL MPI_RECV(ICPACK,(5+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK+1,3,MPI_COMM_WORLD,STATUS,IERROR)
! Unpack the data that we have recieved
        DO K=0,TNKZ
          DO I=0,NKX
            CU1(I,K,NY+1)=ICPACK(I,K,1)
            CU2(I,K,NY+1)=ICPACK(I,K,2)
            CU3(I,K,NY+1)=ICPACK(I,K,3)
            CP(I,K,NY+1)=ICPACK(I,K,4)
            DO N=1,N_TH
              CTH(I,K,NY+1,N)=ICPACK(I,K,4+N)
            END DO
          END DO
        END DO
      END IF

      END IF

      RETURN
      END


      SUBROUTINE GHOST_PHI_MPI
! This subroutine is part of the MPI package for the channel flow
! Diablo package.
! Here, we define a set of ghost cells on each process
! the ghost cells contain information of CR1 from the neighboring nodes
! and allow us to compute finite differences over the local gridpoints.
! We need to update the contents of the ghost cells at the start of
! each Runge-Kutta substep

      include 'header'
      INCLUDE "mpif.h"
      include 'header_mpi'

      integer i,j,k,N

! Define the arrays that will be used for data packing.  This makes the
! communication between processes more efficient by only requiring one
! send and recieve.
! The communication will be done in Fourier space, so these arrays should
! be complex arrays to match the velocity
! The size of the buffer arrys is 0:NKX,0:TNKZ,# of variables
      COMPLEX*16 OCPACK(0:NX/3,0:2*NZ/3,1)
      COMPLEX*16 ICPACK(0:NX/3,0:2*NZ/3,1)

! If we are using more than one processor, then we need to pass data
      IF (NPROCS.gt.1) THEN

! First, Pass data up the chain to higher ranked processes

      IF (RANK.eq.0) THEN
! If we are the lowest ranked process, then we don't need to recieve
! data at the lower ghost cells, these will be filled with boundary
! condition information
        DO K=0,TNKZ
          DO I=0,NKX
            OCPACK(I,K,1)=CR1(I,K,NY)
          END DO
        END DO

!        OCPACK(0,0,5+N_TH)=0.
! Now, we have packed the data into a compact array, pass the data up
        CALL MPI_SEND(OCPACK,(1)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK+1,1,MPI_COMM_WORLD,IERROR)

! End if RANK=0
      ELSE IF (RANK.lt.NPROCS-1) THEN
! Here, we are one of the middle processes and we need to pass data
! up and recieve data from below
        DO K=0,TNKZ
          DO I=0,NKX
            OCPACK(I,K,1)=CR1(I,K,NY)
          END DO
        END DO

!        OCPACK(0,0,5+N_TH)=0.
! Use MPI_SENDRECV since we need to recieve and send data
        CALL MPI_SEND(OCPACK,(1)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK+1,1,MPI_COMM_WORLD,IERROR)

        CALL MPI_RECV(ICPACK,(1)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK-1,1,MPI_COMM_WORLD,STATUS,IERROR)
! Now, unpack the data that we have recieved
        DO K=0,TNKZ
          DO I=0,NKX
            CR1(I,K,1)=ICPACK(I,K,1)
          END DO
        END DO

      ELSE
! Otherwise, we must be the uppermost process with RANK=NPROCS-1
! Here, we need to recieve data from below, but don't need to send data up
        CALL MPI_RECV(ICPACK,(1)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK-1,1,MPI_COMM_WORLD,STATUS,IERROR)
! Unpack the data that we have recieved
        DO K=0,TNKZ
          DO I=0,NKX
            CR1(I,K,1)=ICPACK(I,K,1)
          END DO
        END DO
      END IF

! AT this point we have passed data up the chain

! Now, we need to pass data down the chain of processes
      IF (RANK.eq.NPROCS-1) THEN
! If we are the higest ranked process, then we don't need to recieve
! data at the upper ghost cells, these will be filled with boundary
! condition information
        DO K=0,TNKZ
          DO I=0,NKX
            OCPACK(I,K,1)=CR1(I,K,2)
          END DO
        END DO

!        OCPACK(0,0,5+N_TH)=0.
! Now, we have packed the data into a compact array, pass the data up
        CALL MPI_SEND(OCPACK,(1)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK-1,3,MPI_COMM_WORLD,IERROR)
      ELSE IF (RANK.GT.0) THEN
! Here, we are one of the middle processes and we need to pass data
! down and recieve data from above us
        DO K=0,TNKZ
          DO I=0,NKX
            OCPACK(I,K,1)=CR1(I,K,2)
          END DO
        END DO
!        OCPACK(0,0,5+N_TH)=0.
        CALL MPI_SEND(OCPACK,(1)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK-1,3,MPI_COMM_WORLD,IERROR)

        CALL MPI_RECV(ICPACK,(1)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK+1,3,MPI_COMM_WORLD,STATUS,IERROR)
! Now, unpack the data that we have recieved
        DO K=0,TNKZ
          DO I=0,NKX
            CR1(I,K,NY+1)=ICPACK(I,K,1)
          END DO
        END DO
      ELSE
! Here, we must be the lowest process (RANK=0) and we need to recieve
! data from above
        CALL MPI_RECV(ICPACK,(1)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK+1,3,MPI_COMM_WORLD,STATUS,IERROR)
! Unpack the data that we have recieved
        DO K=0,TNKZ
          DO I=0,NKX
            CR1(I,K,NY+1)=ICPACK(I,K,1)
          END DO
        END DO
      END IF

      END IF

      RETURN
      END




      SUBROUTINE GHOST_ZERO_CHAN_MPI_previous

! This subroutine is part of the MPI package for the channel flow
! Diablo package.
! Here, we define a set of ghost cells on each process
! the ghost cells contain information from the neighboring nodes
! and allow us to compute finite differences over the local gridpoints.
! We need to update the contents of the ghost zero cells at the start of
! les process

      include 'header'
      INCLUDE "mpif.h"
      include 'header_mpi'

      integer i,j,k,N

! Define the arrays that will be used for data packing.  This makes the
! communication between processes more efficient by only requiring one
! send and recieve.
! The communication will be done in Fourier space, so these arrays should
! be complex arrays to match the velocity
! The size of the buffer arrys is 0:NKX,0:TNKZ,# of variables
      COMPLEX*8 OCPACK(0:NX/3,0:2*NZ/3,7+N_TH)
      COMPLEX*8 ICPACK(0:NX/3,0:2*NZ/3,7+N_TH)

! If we are using more than one processor, then we need to pass data
      IF (NPROCS.gt.1) THEN

! First, Pass data up the chain to higher ranked processes

      IF (RANK.eq.0) THEN
! If we are the lowest ranked process, then we don't need to recieve
! data at the lower ghost cells, these will be filled with boundary
! condition information
        DO K=0,TNKZ
          DO I=0,NKX
            OCPACK(I,K,1)=CU1(I,K,NY-1)
            OCPACK(I,K,2)=CU2(I,K,NY-1)
            OCPACK(I,K,3)=CU3(I,K,NY-1)
            OCPACK(I,K,4)=CU1(I,K,NY-2)
            OCPACK(I,K,5)=CU2(I,K,NY-2)
            OCPACK(I,K,6)=CU3(I,K,NY-2)                  
            DO N=1,N_TH
              OCPACK(I,K,6+N)=CTH(I,K,NY-1,N)
              OCPACK(I,K,7+N)=CTH(I,K,NY-2,N)
            END DO
          END DO
        END DO            
! Now, we have packed the data into a compact array, pass the data up
        CALL MPI_SEND(OCPACK,(7+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK+1,1,MPI_COMM_WORLD,IERROR)

! End if RANK=0
      ELSE IF (RANK.lt.NPROCS-1) THEN
! Here, we are one of the middle processes and we need to pass data
! up and recieve data from below
        DO K=0,TNKZ
          DO I=0,NKX
            OCPACK(I,K,1)=CU1(I,K,NY-1)
            OCPACK(I,K,2)=CU2(I,K,NY-1)
            OCPACK(I,K,3)=CU3(I,K,NY-1)
            OCPACK(I,K,4)=CU1(I,K,NY-2)
            OCPACK(I,K,5)=CU2(I,K,NY-2)
            OCPACK(I,K,6)=CU3(I,K,NY-2) 
            DO N=1,N_TH
              OCPACK(I,K,6+N)=CTH(I,K,NY-1,N)
              OCPACK(I,K,7+N)=CTH(I,K,NY-2,N)
            END DO
          END DO
        END DO
! Use MPI_SENDRECV since we need to recieve and send data
        CALL MPI_SEND(OCPACK,(7+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK+1,1,MPI_COMM_WORLD,IERROR)

        CALL MPI_RECV(ICPACK,(7+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK-1,1,MPI_COMM_WORLD,STATUS,IERROR)
! Now, unpack the data that we have recieved
        DO K=0,TNKZ
          DO I=0,NKX
            CU1(I,K,0)=ICPACK(I,K,1)
            CU2(I,K,0)=ICPACK(I,K,2)
            CU3(I,K,0)=ICPACK(I,K,3)
            CU1(I,K,-1)=ICPACK(I,K,4)
            CU2(I,K,-1)=ICPACK(I,K,5)
            CU3(I,K,-1)=ICPACK(I,K,6) 
            DO N=1,N_TH
              CTH(I,K,0,N)=ICPACK(I,K,6+N)
              CTH(I,K,-1,N)=ICPACK(I,K,7+N)
            END DO
          END DO
        END DO

      ELSE
! Otherwise, we must be the uppermost process with RANK=NPROCS-1
! Here, we need to recieve data from below, but don't need to send data up
        CALL MPI_RECV(ICPACK,(7+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK-1,1,MPI_COMM_WORLD,STATUS,IERROR)
! Unpack the data that we have recieved
        DO K=0,TNKZ
          DO I=0,NKX
            CU1(I,K,0)=ICPACK(I,K,1)
            CU2(I,K,0)=ICPACK(I,K,2)
            CU3(I,K,0)=ICPACK(I,K,3)
            CU1(I,K,-1)=ICPACK(I,K,4)
            CU2(I,K,-1)=ICPACK(I,K,5)
            CU3(I,K,-1)=ICPACK(I,K,6)
            DO N=1,N_TH
              CTH(I,K,0,N)=ICPACK(I,K,6+N)
              CTH(I,K,-1,N)=ICPACK(I,K,7+N)
            END DO
          END DO
        END DO
      END IF

      ! Now, we need to pass data down the chain of processes
      IF (RANK.eq.NPROCS-1) THEN
! If we are the higest ranked process, then we don't need to recieve
! data at the upper ghost cells, these will be filled with boundary
! condition information
        DO K=0,TNKZ
          DO I=0,NKX
            OCPACK(I,K,1)=CU1(I,K,3)
            OCPACK(I,K,2)=CU2(I,K,3)
            OCPACK(I,K,3)=CU3(I,K,3)
            DO N=1,N_TH
              OCPACK(I,K,3+N)=CTH(I,K,3,N)
            END DO
          END DO
        END DO
! Now, we have packed the data into a compact array, pass the data up
        CALL MPI_SEND(OCPACK,(3+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK-1,3,MPI_COMM_WORLD,IERROR)
      ELSE IF (RANK.GT.0) THEN
! Here, we are one of the middle processes and we need to pass data
! down and recieve data from above us
        DO K=0,TNKZ
          DO I=0,NKX
            OCPACK(I,K,1)=CU1(I,K,3)
            OCPACK(I,K,2)=CU2(I,K,3)
            OCPACK(I,K,3)=CU3(I,K,3)
            DO N=1,N_TH
              OCPACK(I,K,3+N)=CTH(I,K,3,N)
            END DO
          END DO
        END DO

        CALL MPI_SEND(OCPACK,(3+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK-1,3,MPI_COMM_WORLD,IERROR)

        CALL MPI_RECV(ICPACK,(3+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK+1,3,MPI_COMM_WORLD,STATUS,IERROR)
! Now, unpack the data that we have recieved
        DO K=0,TNKZ
          DO I=0,NKX
            CU1(I,K,NY+2)=ICPACK(I,K,1)
            CU2(I,K,NY+2)=ICPACK(I,K,2)
            CU3(I,K,NY+2)=ICPACK(I,K,3)
            DO N=1,N_TH
              CTH(I,K,NY+2,N)=ICPACK(I,K,3+N)
            END DO
          END DO
        END DO
      ELSE
! Here, we must be the lowest process (RANK=0) and we need to recieve
! data from above
        CALL MPI_RECV(ICPACK,(3+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK+1,3,MPI_COMM_WORLD,STATUS,IERROR)
! Unpack the data that we have recieved
        DO K=0,TNKZ
          DO I=0,NKX
            CU1(I,K,NY+2)=ICPACK(I,K,1)
            CU2(I,K,NY+2)=ICPACK(I,K,2)
            CU3(I,K,NY+2)=ICPACK(I,K,3)
            DO N=1,N_TH
              CTH(I,K,NY+2,N)=ICPACK(I,K,3+N)
            END DO
          END DO
        END DO
      END IF

! AT this point we have passed data up the chain

      END IF

      RETURN
      END

      SUBROUTINE GHOST_ZERO_CHAN_MPI

! This subroutine is part of the MPI package for the channel flow
! Diablo package.
! Here, we define a set of ghost cells on each process
! the ghost cells contain information from the neighboring nodes
! and allow us to compute finite differences over the local gridpoints.
! We need to update the contents of the ghost zero cells at the start of
! les process

      include 'header'
      INCLUDE "mpif.h"
      include 'header_mpi'

      integer i,j,k,N

! Define the arrays that will be used for data packing.  This makes the
! communication between processes more efficient by only requiring one
! send and recieve.
! The communication will be done in Fourier space, so these arrays should
! be complex arrays to match the velocity
! The size of the buffer arrys is 0:NKX,0:TNKZ,# of variables
      COMPLEX*16 OCPACK(0:NX/3,0:2*NZ/3,7+N_TH+1)
      COMPLEX*16 ICPACK(0:NX/3,0:2*NZ/3,7+N_TH+1)

! If we are using more than one processor, then we need to pass data
      IF (NPROCS.gt.1) THEN

! First, Pass data up the chain to higher ranked processes

      IF (RANK.eq.0) THEN
! If we are the lowest ranked process, then we don't need to recieve
! data at the lower ghost cells, these will be filled with boundary
! condition information
        DO K=0,TNKZ
          DO I=0,NKX
            OCPACK(I,K,1)=CU1(I,K,NY-1)
            OCPACK(I,K,2)=CU2(I,K,NY-1)
            OCPACK(I,K,3)=CU3(I,K,NY-1)
            OCPACK(I,K,4)=CU1(I,K,NY-2)
            OCPACK(I,K,5)=CU2(I,K,NY-2)
            OCPACK(I,K,6)=CU3(I,K,NY-2)                  
            DO N=1,N_TH
              OCPACK(I,K,6+N)=CTH(I,K,NY-1,N)
              OCPACK(I,K,7+N)=CTH(I,K,NY-2,N)
            END DO
          END DO
        END DO
        OCPACK(0,0,8+N_TH) = 0.           
! Now, we have packed the data into a compact array, pass the data up
        CALL MPI_SEND(OCPACK,(8+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK+1,1,MPI_COMM_WORLD,IERROR)

! End if RANK=0
      ELSE IF (RANK.lt.NPROCS-1) THEN
! Here, we are one of the middle processes and we need to pass data
! up and recieve data from below
        DO K=0,TNKZ
          DO I=0,NKX
            OCPACK(I,K,1)=CU1(I,K,NY-1)
            OCPACK(I,K,2)=CU2(I,K,NY-1)
            OCPACK(I,K,3)=CU3(I,K,NY-1)
            OCPACK(I,K,4)=CU1(I,K,NY-2)
            OCPACK(I,K,5)=CU2(I,K,NY-2)
            OCPACK(I,K,6)=CU3(I,K,NY-2) 
            DO N=1,N_TH
              OCPACK(I,K,6+N)=CTH(I,K,NY-1,N)
              OCPACK(I,K,7+N)=CTH(I,K,NY-2,N)
            END DO
          END DO
        END DO
        OCPACK(0,0,8+N_TH) = 0.
! Use MPI_SENDRECV since we need to recieve and send data
        CALL MPI_SEND(OCPACK,(8+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK+1,1,MPI_COMM_WORLD,IERROR)

        CALL MPI_RECV(ICPACK,(8+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK-1,1,MPI_COMM_WORLD,STATUS,IERROR)
! Now, unpack the data that we have recieved
        DO K=0,TNKZ
          DO I=0,NKX
            CU1(I,K,0)=ICPACK(I,K,1)
            CU2(I,K,0)=ICPACK(I,K,2)
            CU3(I,K,0)=ICPACK(I,K,3)
            CU1(I,K,-1)=ICPACK(I,K,4)
            CU2(I,K,-1)=ICPACK(I,K,5)
            CU3(I,K,-1)=ICPACK(I,K,6) 
            DO N=1,N_TH
              CTH(I,K,0,N)=ICPACK(I,K,6+N)
              CTH(I,K,-1,N)=ICPACK(I,K,7+N)
            END DO
          END DO
        END DO

      ELSE
! Otherwise, we must be the uppermost process with RANK=NPROCS-1
! Here, we need to recieve data from below, but don't need to send data up
        CALL MPI_RECV(ICPACK,(8+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK-1,1,MPI_COMM_WORLD,STATUS,IERROR)
! Unpack the data that we have recieved
        DO K=0,TNKZ
          DO I=0,NKX
            CU1(I,K,0)=ICPACK(I,K,1)
            CU2(I,K,0)=ICPACK(I,K,2)
            CU3(I,K,0)=ICPACK(I,K,3)
            CU1(I,K,-1)=ICPACK(I,K,4)
            CU2(I,K,-1)=ICPACK(I,K,5)
            CU3(I,K,-1)=ICPACK(I,K,6)
            DO N=1,N_TH
              CTH(I,K,0,N)=ICPACK(I,K,6+N)
              CTH(I,K,-1,N)=ICPACK(I,K,7+N)
            END DO
          END DO
        END DO
      END IF

      ! Now, we need to pass data down the chain of processes
      IF (RANK.eq.NPROCS-1) THEN
! If we are the higest ranked process, then we don't need to recieve
! data at the upper ghost cells, these will be filled with boundary
! condition information
        DO K=0,TNKZ
          DO I=0,NKX
            OCPACK(I,K,1)=CU1(I,K,3)
            OCPACK(I,K,2)=CU2(I,K,3)
            OCPACK(I,K,3)=CU3(I,K,3)
            DO N=1,N_TH
              OCPACK(I,K,3+N)=CTH(I,K,3,N)
            END DO
          END DO
        END DO
        OCPACK(0,0,4+N_TH) = 0.
! Now, we have packed the data into a compact array, pass the data up
        CALL MPI_SEND(OCPACK,(4+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK-1,3,MPI_COMM_WORLD,IERROR)
      ELSE IF (RANK.GT.0) THEN
! Here, we are one of the middle processes and we need to pass data
! down and recieve data from above us
        DO K=0,TNKZ
          DO I=0,NKX
            OCPACK(I,K,1)=CU1(I,K,3)
            OCPACK(I,K,2)=CU2(I,K,3)
            OCPACK(I,K,3)=CU3(I,K,3)
            DO N=1,N_TH
              OCPACK(I,K,3+N)=CTH(I,K,3,N)
            END DO
          END DO
        END DO
        OCPACK(0,0,4+N_TH) = 0.
        CALL MPI_SEND(OCPACK,(4+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK-1,3,MPI_COMM_WORLD,IERROR)

        CALL MPI_RECV(ICPACK,(4+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK+1,3,MPI_COMM_WORLD,STATUS,IERROR)
! Now, unpack the data that we have recieved
        DO K=0,TNKZ
          DO I=0,NKX
            CU1(I,K,NY+2)=ICPACK(I,K,1)
            CU2(I,K,NY+2)=ICPACK(I,K,2)
            CU3(I,K,NY+2)=ICPACK(I,K,3)
            DO N=1,N_TH
              CTH(I,K,NY+2,N)=ICPACK(I,K,3+N)
            END DO
          END DO
        END DO
      ELSE
! Here, we must be the lowest process (RANK=0) and we need to recieve
! data from above
        CALL MPI_RECV(ICPACK,(4+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK+1,3,MPI_COMM_WORLD,STATUS,IERROR)
! Unpack the data that we have recieved
        DO K=0,TNKZ
          DO I=0,NKX
            CU1(I,K,NY+2)=ICPACK(I,K,1)
            CU2(I,K,NY+2)=ICPACK(I,K,2)
            CU3(I,K,NY+2)=ICPACK(I,K,3)
            DO N=1,N_TH
              CTH(I,K,NY+2,N)=ICPACK(I,K,3+N)
            END DO
          END DO
        END DO
      END IF

! AT this point we have passed data up the chain

      END IF

      RETURN
      END

      SUBROUTINE GHOST_GRID_MPI
! This subroutine is part of the MPI package for the channel flow
! Diablo package.
! Here, we define a set of ghost cells on each process
! the ghost cells contain information from the neighboring nodes
! and allow us to compute finite differences over the local gridpoints.
! We need to update the contents of the ghost cells at the start of
! each Runge-Kutta substep

      include 'header'
      INCLUDE "mpif.h"

      include 'header_mpi'

      integer i,j,k,N

      real*8 OCPACK(5),ICPACK(5)

! First, Pass data up the chain to higher ranked processes

      IF (RANK.eq.0) THEN

! Set the lower ghost cells
        GYF(0) =2.d0*GYF(1)-GYF(2)
        GY(0)  =2.d0*GY(1)-GY(2)
        GYF(-1)=2.d0*GYF(0)-GYF(1)

        OCPACK(1)=GYF(NY-1)
        OCPACK(2)=GYF(NY)
        OCPACK(3)=GY(NY)       
        OCPACK(4)=GY(NY-1)
        OCPACK(5)=GYF(NY-2)

        CALL MPI_SEND(OCPACK,5
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANK+1,1,MPI_COMM_WORLD,IERROR)

! End if RANK=0
      ELSE IF (RANK.lt.NPROCS-1) THEN
! Here, we are one of the middle processes and we need to pass data
! up and recieve data from below
        OCPACK(1)=GYF(NY-1)
        OCPACK(2)=GYF(NY)
        OCPACK(3)=GY(NY)
        OCPACK(4)=GY(NY-1)
        OCPACK(5)=GYF(NY-2)

        CALL MPI_SEND(OCPACK,5
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANK+1,1,MPI_COMM_WORLD,IERROR)

        CALL MPI_RECV(ICPACK,5
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANK-1,1,MPI_COMM_WORLD,STATUS,IERROR)
! Now, unpack the data that we have recieved
        GYF(0) =ICPACK(1)
        GYF(1) =ICPACK(2)
        GY(1)  =ICPACK(3)
        GY(0)  =ICPACK(4)
        GYF(-1)=ICPACK(5)

      ELSE
! Otherwise, we must be the uppermost process with RANK=NPROCS-1
! Set the top ghost cell
        GYF(NY+1)=2.*GYF(NY)-GYF(NYM)
        GYF(NY+2)=2.*GYF(NY+1)-GYF(NY)
        GY(NY+2) =2.*GY(NY+1)-GY(NY)

! Here, we need to recieve data from below, but don't need to send data up
        CALL MPI_RECV(ICPACK,5
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANK-1,1,MPI_COMM_WORLD,STATUS,IERROR)
! Now, unpack the data that we have recieved
        GYF(0) =ICPACK(1)
        GYF(1) =ICPACK(2)
        GY(1)  =ICPACK(3)
        GY(0)  =ICPACK(4)
        GYF(-1)=ICPACK(5)

      END IF
! AT this point we have passed data up the chain

      IF (RANK.eq.NPROCS-1) THEN
! If we are the higest ranked process, then we don't need to recieve
! data at the upper ghost cells, these will be filled with boundary
! condition information

        OCPACK(1)=GYF(2)
        OCPACK(2)=GY(2)
        OCPACK(3)=GYF(3)
        OCPACK(4)=GY(3)

! Now, we have packed the data into a compact array, pass the data up
        CALL MPI_SEND(OCPACK,4
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANK-1,3,MPI_COMM_WORLD,IERROR)
      ELSE IF (RANK.GT.0) THEN
! Here, we are one of the middle processes and we need to pass data
! down and recieve data from above us

        OCPACK(1)=GYF(2)
        OCPACK(2)=GY(2)
        OCPACK(3)=GYF(3)
        OCPACK(4)=GY(3)
        CALL MPI_SEND(OCPACK,4
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANK-1,3,MPI_COMM_WORLD,IERROR)
        CALL MPI_RECV(ICPACK,4
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANK+1,3,MPI_COMM_WORLD,STATUS,IERROR)
! Now, unpack the data that we have recieved
        GYF(NY+1)=ICPACK(1)
        GY(NY+1) =ICPACK(2)
        GYF(NY+2)=ICPACK(3)
        GY(NY+2) =ICPACK(4)
      ELSE
! Here, we must be the lowest process (RANK=0) and we need to recieve
! data from above
        CALL MPI_RECV(ICPACK,4
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANK+1,3,MPI_COMM_WORLD,STATUS,IERROR)
! Unpack the data that we have recieved
        GYF(NY+1)=ICPACK(1)
        GY(NY+1) =ICPACK(2)
        GYF(NY+2)=ICPACK(3)
        GY(NY+2) =ICPACK(4)
      END IF

      RETURN
      END


      SUBROUTINE MPI_BCAST_REAL(IN_R,Z_SIZE,Y_SIZE)
 
      include 'header'
      INCLUDE "mpif.h"

      include 'header_mpi'

      integer i,j,k,N,send_rank,Z_SIZE,Y_SIZE
      real*8 IN_R(0:Z_SIZE*Y_SIZE-1) 

      send_rank = 0


      CALL MPI_BCAST(IN_R,Z_SIZE*Y_SIZE,MPI_DOUBLE_PRECISION,
     &send_rank,MPI_COMM_WORLD,IERROR)

      RETURN
      END

      SUBROUTINE MPI_BCAST_REAL_L(IN_R,Z_SIZE,Y_SIZE)
 
      include 'header'
      INCLUDE "mpif.h"

      include 'header_mpi'

      integer i,j,k,N,send_rank,Z_SIZE,Y_SIZE
      real*8 IN_R(0:Z_SIZE*Y_SIZE-1) 

      send_rank = NPROCS-1


      CALL MPI_BCAST(IN_R,Z_SIZE*Y_SIZE,MPI_DOUBLE_PRECISION,
     &send_rank,MPI_COMM_WORLD,IERROR)

      RETURN
      END
 
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE MPI_SECLT_MAX(OPACK)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      include 'header'
      INCLUDE "mpif.h"

      include 'header_mpi'


       REAL(8),DIMENSION(1:NPROCS) :: IPACK
       REAL(8) :: OPACK

! READ IN LOCAL DT BASED ON CFL
!      OPACK = DT

! SEND ALL LOCAL DT'S TO EVERY PROCESS
       CALL MPI_ALLGATHER(OPACK,1,MPI_DOUBLE_PRECISION,IPACK,1,
     & MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERROR)

       OPACK = MAXVAL(IPACK)

       RETURN
       END


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE MPI_SECLT(OPACK)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      include 'header'
      INCLUDE "mpif.h"

      include 'header_mpi'


       REAL(8),DIMENSION(1:NPROCS) :: IPACK
       REAL(8) :: OPACK

! READ IN LOCAL DT BASED ON CFL
!      OPACK = DT

! SEND ALL LOCAL DT'S TO EVERY PROCESS
       CALL MPI_ALLGATHER(OPACK,1,MPI_DOUBLE_PRECISION,IPACK,1,
     & MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERROR)

       OPACK = MINVAL(IPACK)

       RETURN
       END


      SUBROUTINE COURANT_MPI
! This subroutine is part of the MPI package for the channel flow
! Diablo package.

      include 'header'
      INCLUDE "mpif.h"

      include 'header_mpi'

      integer i,j,k,N

      real*8 OCPACK,ICPACK,dt_temp

! First, Pass data up the chain to higher ranked processes

      IF (RANK.eq.0) THEN

! Set the dt_temp according to the lower most processor
! send it to the higher ranked processor for comparison

        dt_temp = dt
        OCPACK=dt_temp

        CALL MPI_SEND(OCPACK,1
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANK+1,1,MPI_COMM_WORLD,IERROR)

! End if RANK=0
      ELSE IF (RANK.lt.NPROCS-1) THEN
! Here, we are one of the middle processes and we need to 
! recieve data from below and pass data up  
       
        CALL MPI_RECV(ICPACK,1
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANK-1,1,MPI_COMM_WORLD,STATUS,IERROR)
! Now, unpack the data that we have recieved
        dt_temp=ICPACK

        IF (dt .le. dt_temp) THEN
          dt_temp = dt
        ENDIF
        
        OCPACK=dt_temp        

        CALL MPI_SEND(OCPACK,1
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANK+1,1,MPI_COMM_WORLD,IERROR)


      ELSE
! Otherwise, we must be the uppermost process with RANK=NPROCS-1


! Here, we need to recieve data from below, but don't need to send data up
        CALL MPI_RECV(ICPACK,1
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANK-1,1,MPI_COMM_WORLD,STATUS,IERROR)
! Now, unpack the data that we have recieved

        dt_temp=ICPACK

        IF (dt .le. dt_temp) THEN
          dt_temp = dt
        ENDIF
      END IF
! AT this point we have passed data up the chain

      IF (RANK.eq.NPROCS-1) THEN
! If we are the higest ranked process, then we don't need to recieve

        dt = dt_temp
        OCPACK=dt_temp
! Now, we have packed the data into a compact array, pass the data up
        CALL MPI_SEND(OCPACK,1
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANK-1,3,MPI_COMM_WORLD,IERROR)
      ELSE 
! Here, we are one of the middle processes and we need to pass data
! down and recieve data from above us

        CALL MPI_RECV(ICPACK,1
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANK+1,3,MPI_COMM_WORLD,STATUS,IERROR)
! Now, unpack the data that we have recieved
        dt=ICPACK

      END IF

      RETURN
      END

      SUBROUTINE COURANT_MPI_1
! This subroutine is part of the MPI package for the channel flow
! Diablo package.

      include 'header'
      INCLUDE "mpif.h"

      include 'header_mpi'

      integer i,j,k,N

      real*8 OCPACK,ICPACK,dt_temp

     
    
      IF (RANK.NE.0) THEN
! Here, we are one of the middle processes and we need to 
! recieve data from below and pass data up  
       
        CALL MPI_RECV(ICPACK,1
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANK-1,1,MPI_COMM_WORLD,STATUS,IERROR)
! Now, unpack the data that we have recieved
        dt_temp=ICPACK

        IF (dt .le. dt_temp) THEN
          dt_temp = dt
        ENDIF
        
        IF (RANK.NE.NPROCS-1) THEN

          OCPACK=dt_temp        

          CALL MPI_SEND(OCPACK,1
     &                 ,MPI_DOUBLE_PRECISION
     &                 ,RANK+1,1,MPI_COMM_WORLD,IERROR)
        ELSE
          dt = dt_temp
          OCPACK=dt_temp

        DO I = 1,NPROCS-1
         CALL MPI_SEND(OCPACK,1
     &               ,MPI_DOUBLE_PRECISION
     &               ,I-1,1,MPI_COMM_WORLD,IERROR)
        ENDDO

        ENDIF


      ELSE
        dt_temp = dt
        OCPACK=dt_temp

        CALL MPI_SEND(OCPACK,1
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANK+1,1,MPI_COMM_WORLD,IERROR)
     
      ENDIF


! AT this point we have passed data up the chain

      IF (RANK.NE.NPROCS-1) THEN
 
! Here, we are one of the middle processes and we need to pass data
! down and recieve data from above us

        CALL MPI_RECV(ICPACK,1
     &               ,MPI_DOUBLE_PRECISION
     &               ,NPROCS-1,1,MPI_COMM_WORLD,STATUS,IERROR)
! Now, unpack the data that we have recieved
        dt=ICPACK

      END IF

      RETURN
      END

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE MPI_INST_PROF(OPACK,IPACK)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      include 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'

      REAL(8),DIMENSION(1:NY*NPROCS) :: IPACK
      REAL(8),DIMENSION(1:NY) :: OPACK

! READ IN LOCAL DT BASED ON CFL
!      OPACK = DT

! SEND ALL LOCAL DT'S TO EVERY PROCESS
      CALL MPI_ALLGATHER(OPACK,NY,MPI_DOUBLE_PRECISION,IPACK,NY, 
     + MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERROR)

!      DT = MINVAL(IPACK)

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE THOMAS_FORWARD_REAL_MPI(A,B,C,G,NY,NX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This version solves for one full row, then passes the data
C This subroutine performs the backward sweep of the Thomas algorithm
C Thomas algorithm solves Ax=b for tridiagonal A
C The RHS vector and solution are real
C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
C Returns solution in x
C The indexing should be done by ROW, ie.
C [ b1  c1   0   0   0 ...
C [ a2  b2  c2   0   0 ...
C [  0  a3  b3   c3  0 ...

      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'

      INTEGER I,J,NX,NY
      REAL*8 A(0:NX-1,0:NY+1), B(0:NX-1,0:NY+1), C(0:NX-1,0:NY+1)
     &         , G(0:NX-1,0:NY+1)

      REAL*8 OCPACK(0:NX-1,4),ICPACK(0:NX-1,4)


      IF (RANK.NE.0) THEN
C If we aren't the lowest process, then wait for data
        CALL MPI_RECV(OCPACK,4*NX,MPI_DOUBLE_PRECISION,RANK-1,12
     &               ,MPI_COMM_WORLD,status,ierror)
C Unpack the data
        DO I=0,NX-1
        A(I,1)=OCPACK(I,1)
        B(I,1)=OCPACK(I,2)
        C(I,1)=OCPACK(I,3)
        G(I,1)=OCPACK(I,4)
        END DO
      ELSE
! We are the lowest process, so solve at J=1
        J=1
        DO I=0,NX-1
        A(I,J)=-A(I,J)/B(I,J-1)
        B(I,J)=B(I,J)+A(I,J)*C(I,J-1)
        G(I,J)=G(I,J)+A(I,J)*G(I,J-1)
        ENDDO
       END IF

! For all processes, solve from 2..NY 
        DO J=2,NY+1
          DO I=0,NX-1
            A(I,J)=-A(I,J)/B(I,J-1)
            B(I,J)=B(I,J)+A(I,J)*C(I,J-1)
            G(I,J)=G(I,J)+A(I,J)*G(I,J-1)
          END DO
        END DO

      IF (RANK.NE.NPROCS-1) THEN
        DO I=0,NX-1
        ICPACK(I,1)=A(I,NY)
        ICPACK(I,2)=B(I,NY)
        ICPACK(I,3)=C(I,NY)
        ICPACK(I,4)=G(I,NY)
        END DO
        CALL MPI_SEND(ICPACK,4*NX,MPI_DOUBLE_PRECISION,RANK+1,12
     &               ,MPI_COMM_WORLD,ierror)
      END IF

      RETURN
      END



C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE THOMAS_FORWARD_COMPLEX_MPI(A,B,C,G,NY,NX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
C This subroutine performs the backward sweep of the Thomas algorithm
C Thomas algorithm solves Ax=b for tridiagonal A
C The RHS vector and solution are real
C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
C Returns solution in x
C The indexing should be done by ROW, ie.
C [ b1  c1   0   0   0 ...
C [ a2  b2  c2   0   0 ...
C [  0  a3  b3   c3  0 ...

      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'

      INTEGER I,J,NX,NY
      REAL*8 A(0:NX,0:NY+1), B(0:NX,0:NY+1), C(0:NX,0:NY+1)
      COMPLEX*16  G(0:NX,0:NY+1)

      COMPLEX*16 OCPACK(0:NX,4),ICPACK(0:NX,4)

c      i=10
c      do j=0,NY+1
c        write(700+rank,*) j,A(i,j),B(i,j),C(i,j),G(i,j)
c      end do

      IF (RANK.NE.0) THEN
C If we aren't the lowest process, then wait for data
        CALL MPI_RECV(OCPACK,4*(NX+1),MPI_DOUBLE_COMPLEX,RANK-1,13
     &               ,MPI_COMM_WORLD,status,ierror)
C Unpack the data
        DO I=0,NX
        A(I,1)=real(OCPACK(I,1))
        B(I,1)=real(OCPACK(I,2))
        C(I,1)=real(OCPACK(I,3))
        G(I,1)=OCPACK(I,4)
        END DO
       END IF

! For all processes, solve from 2..NY 
        DO J=2,NY
          DO I=0,NX
            A(I,J)=-A(I,J)/B(I,J-1)
            B(I,J)=B(I,J)+A(I,J)*C(I,J-1)
            G(I,J)=G(I,J)+A(I,J)*G(I,J-1)
          END DO
        END DO

      IF (RANK.NE.NPROCS-1) THEN
        DO I=0,NX
        ICPACK(I,1)=A(I,NY)
        ICPACK(I,2)=B(I,NY)
        ICPACK(I,3)=C(I,NY)
        ICPACK(I,4)=G(I,NY)
        END DO
        CALL MPI_SEND(ICPACK,4*(NX+1),MPI_DOUBLE_COMPLEX,RANK+1,13
     &               ,MPI_COMM_WORLD,ierror)
      END IF

      RETURN
      END




C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE THOMAS_BACKWARD_REAL_MPI(A,B,C,G,NY,NX)
C----*|--.---------.---------.---------.---------.---------.---------.-|------
C This subroutine performs the backward sweep of the Thomas algorithm
C Thomas algorithm solves Ax=b for tridiagonal A
C The RHS vector and solution are real
C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
C Returns solution in x
C The indexing should be done by ROW, ie.
C [ b1  c1   0   0   0 ...
C [ a2  b2  c2   0   0 ...
C [  0  a3  b3   c3  0 ...

      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'

      INTEGER I,J,NX,NY
      REAL*8 A(0:NX-1,0:NY+1), B(0:NX-1,0:NY+1), C(0:NX-1,0:NY+1)
     &     , G(0:NX-1,0:NY+1)
      REAL*8 ICPACK(0:NX-1),OCPACK(0:NX-1)

      IF (RANK.NE.NPROCS-1) THEN
C If we aren't the highest process, then wait for data
        CALL MPI_RECV(OCPACK,NX,MPI_DOUBLE_PRECISION,RANK+1,10
     &               ,MPI_COMM_WORLD,status,ierror)
        DO I=0,NX-1
          G(I,NY+1)=OCPACK(I)
        END DO
      ELSE
C Else, if we are the highest process, compute the solution at j=NY
       DO I=0,NX-1  
         G(I,NY+1)=G(I,NY+1)/B(I,NY+1)
       END DO
      END IF

C All processes solve from NY..1
      DO J=NY,0,-1
      DO I=0,NX-1
        G(I,J)=(G(I,J)-C(I,J)*G(I,J+1))/B(I,J)
      END DO
      END DO

      IF (RANK.NE.0) THEN
        DO I=0,NX-1
          ICPACK(I)=G(I,2)
        END DO
        CALL MPI_SEND(ICPACK,NX,MPI_DOUBLE_PRECISION,RANK-1,10
     &               ,MPI_COMM_WORLD,ierror)
      END IF

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE THOMAS_BACKWARD_COMPLEX_MPI(A,B,C,G,NY,NX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
C This subroutine performs the backward sweep of the Thomas algorithm
C Thomas algorithm solves Ax=b for tridiagonal A
C The RHS vector and solution are real
C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
C Returns solution in x
C The indexing should be done by ROW, ie.
C [ b1  c1   0   0   0 ...
C [ a2  b2  c2   0   0 ...
C [  0  a3  b3   c3  0 ...

      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'

      INTEGER I,J,NX,NY
      REAL*8 A(0:NX,0:NY+1), B(0:NX,0:NY+1), C(0:NX,0:NY+1)
      COMPLEX*16 G(0:NX,0:NY+1)
      COMPLEX*16 ICPACK(0:NX),OCPACK(0:NX)

      IF (RANK.NE.NPROCS-1) THEN
C If we aren't the highest process, then wait for data
        CALL MPI_RECV(OCPACK,NX+1,MPI_DOUBLE_COMPLEX,RANK+1,11
     &               ,MPI_COMM_WORLD,status,ierror)
        DO I=0,NX
          G(I,NY+1)=OCPACK(I)
        END DO
C Solve as usual at j=NY
        J=NY
        do I=0,NX
          G(I,J)=(G(I,J)-C(I,J)*G(I,J+1))/B(I,J)
        end do
      ELSE
C Else, if we are the highest process, compute the solution directly at j=NY
        DO I=0,NX
          G(I,NY)=G(I,NY)/B(I,NY)
        END DO
      END IF

C All processes solve from NY-1..1
      DO J=NY-1,1,-1
      DO I=0,NX
        G(I,J)=(G(I,J)-C(I,J)*G(I,J+1))/B(I,J)
      END DO
      END DO

      IF (RANK.NE.0) THEN
        DO I=0,NX
          ICPACK(I)=G(I,2)
        END DO
        CALL MPI_SEND(ICPACK,NX+1,MPI_DOUBLE_COMPLEX,RANK-1,11
     &               ,MPI_COMM_WORLD,ierror)
      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE THOMAS_FORWARD_REAL_MPI_old(A,B,C,G,NY,NX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This version solves for one full row, then passes the data
C This subroutine performs the backward sweep of the Thomas algorithm
C Thomas algorithm solves Ax=b for tridiagonal A
C The RHS vector and solution are real
C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
C Returns solution in x
C The indexing should be done by ROW, ie.
C [ b1  c1   0   0   0 ...
C [ a2  b2  c2   0   0 ...
C [  0  a3  b3   c3  0 ...

      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'

      INTEGER I,J,NX,NY
      REAL*8 A(0:NX-1,0:NY+1), B(0:NX-1,0:NY+1), C(0:NX-1,0:NY+1)
     &         , G(0:NX-1,0:NY+1)

      REAL*8 OCPACK(0:NX-1,4),ICPACK(0:NX-1,4)


      IF (RANK.NE.0) THEN
C If we aren't the lowest process, then wait for data
        CALL MPI_RECV(OCPACK,4*NX,MPI_DOUBLE_PRECISION,RANK-1,12
     &               ,MPI_COMM_WORLD,status,ierror)
C Unpack the data
        DO I=0,NX-1
        A(I,1)=OCPACK(I,1)
        B(I,1)=OCPACK(I,2)
        C(I,1)=OCPACK(I,3)
        G(I,1)=OCPACK(I,4)
        END DO
      ELSE
! We are the lowest process, so solve at J=1
        J=1
        DO I=0,NX-1
          A(I,J)=-A(I,J)/B(I,J-1)
          B(I,J)=B(I,J)+A(I,J)*C(I,J-1)
          G(I,J)=G(I,J)+A(I,J)*G(I,J-1)
        ENDDO
       END IF

! For all processes, solve from 2..NY 
        DO J=2,NY+1
          DO I=0,NX-1
            A(I,J)=-A(I,J)/B(I,J-1)
            B(I,J)=B(I,J)+A(I,J)*C(I,J-1)
            G(I,J)=G(I,J)+A(I,J)*G(I,J-1)
          END DO
        END DO

      IF (RANK.NE.NPROCS-1) THEN
        DO I=0,NX-1
        ICPACK(I,1)=A(I,NY)
        ICPACK(I,2)=B(I,NY)
        ICPACK(I,3)=C(I,NY)
        ICPACK(I,4)=G(I,NY)
        END DO
        CALL MPI_SEND(ICPACK,4*NX,MPI_DOUBLE_PRECISION,RANK+1,12
     &               ,MPI_COMM_WORLD,ierror)
      END IF

      RETURN
      END



C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE THOMAS_FORWARD_COMPLEX_MPI_old(A,B,C,G,NY,NX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
C This subroutine performs the backward sweep of the Thomas algorithm
C Thomas algorithm solves Ax=b for tridiagonal A
C The RHS vector and solution are real
C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
C Returns solution in x
C The indexing should be done by ROW, ie.
C [ b1  c1   0   0   0 ...
C [ a2  b2  c2   0   0 ...
C [  0  a3  b3   c3  0 ...

      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'

      INTEGER I,J,NX,NY
      REAL*8 A(0:NX,0:NY+1), B(0:NX,0:NY+1), C(0:NX,0:NY+1)
      COMPLEX*16  G(0:NX,0:NY+1)

      COMPLEX*16 OCPACK(0:NX,4),ICPACK(0:NX,4)


      IF (RANK.NE.0) THEN
C If we aren't the lowest process, then wait for data
        CALL MPI_RECV(OCPACK,4*(NX+1),MPI_DOUBLE_COMPLEX,RANK-1,13
     &               ,MPI_COMM_WORLD,status,ierror)
C Unpack the data
        DO I=0,NX
        A(I,1)=real(OCPACK(I,1))
        B(I,1)=real(OCPACK(I,2))
        C(I,1)=real(OCPACK(I,3))
        G(I,1)=OCPACK(I,4)
        END DO
       ELSE
! We are the lowest process, so solve at J=0
        J=1
        DO I=0,NX
        A(I,J)=-A(I,J)/B(I,J-1)
        B(I,J)=B(I,J)+A(I,J)*C(I,J-1)
        G(I,J)=G(I,J)+A(I,J)*G(I,J-1)
        END DO
       END IF

! For all processes, solve from 2..NY 
        DO J=2,NY+1
          DO I=0,NX
            A(I,J)=-A(I,J)/B(I,J-1)
            B(I,J)=B(I,J)+A(I,J+1)*C(I,J-1)
            G(I,J)=G(I,J)+A(I,J+1)*G(I,J-1)
          END DO
        END DO

      IF (RANK.NE.NPROCS-1) THEN
        DO I=0,NX
        ICPACK(I,1)=A(I,NY)
        ICPACK(I,2)=B(I,NY)
        ICPACK(I,3)=C(I,NY)
        ICPACK(I,4)=G(I,NY)
        END DO
        CALL MPI_SEND(ICPACK,4*(NX+1),MPI_DOUBLE_COMPLEX,RANK+1,13
     &               ,MPI_COMM_WORLD,ierror)
      END IF

      RETURN
      END




C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE THOMAS_BACKWARD_REAL_MPI_old(A,B,C,G,NY,NX)
C----*|--.---------.---------.---------.---------.---------.---------.-|------
C This subroutine performs the backward sweep of the Thomas algorithm
C Thomas algorithm solves Ax=b for tridiagonal A
C The RHS vector and solution are real
C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
C Returns solution in x
C The indexing should be done by ROW, ie.
C [ b1  c1   0   0   0 ...
C [ a2  b2  c2   0   0 ...
C [  0  a3  b3   c3  0 ...

      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'

      INTEGER I,J,NX,NY
      REAL*8 A(0:NX-1,0:NY+1), B(0:NX-1,0:NY+1), C(0:NX-1,0:NY+1)
     &     , G(0:NX-1,0:NY+1)
      REAL*8 ICPACK(0:NX-1),OCPACK(0:NX-1)

      IF (RANK.NE.NPROCS-1) THEN
C If we aren't the highest process, then wait for data
        CALL MPI_RECV(OCPACK,NX,MPI_DOUBLE_PRECISION,RANK+1,10
     &               ,MPI_COMM_WORLD,status,ierror)
        DO I=0,NX-1
          G(I,NY+1)=OCPACK(I)
        END DO
      ELSE
C Else, if we are the highest process, compute the solution at j=NY
        DO I=0,NX-1
          G(I,NY+1)=G(I,NY+1)/B(I,NY+1)
        END DO
      END IF

C All processes solve from NY..1
      DO J=NY,0,-1
      DO I=0,NX-1
        G(I,J)=(G(I,J)-C(I,J)*G(I,J+1))/B(I,J)
      END DO
      END DO

      IF (RANK.NE.0) THEN
        DO I=0,NX-1
          ICPACK(I)=G(I,2)
        END DO
        CALL MPI_SEND(ICPACK,NX,MPI_DOUBLE_PRECISION,RANK-1,10
     &               ,MPI_COMM_WORLD,ierror)
      END IF

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE THOMAS_BACKWARD_COMPLEX_MPI_old(A,B,C,G,NY,NX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
C This subroutine performs the backward sweep of the Thomas algorithm
C Thomas algorithm solves Ax=b for tridiagonal A
C The RHS vector and solution are real
C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
C Returns solution in x
C The indexing should be done by ROW, ie.
C [ b1  c1   0   0   0 ...
C [ a2  b2  c2   0   0 ...
C [  0  a3  b3   c3  0 ...

      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'

      INTEGER I,J,NX,NY
      REAL*8 A(0:NX,0:NY+1), B(0:NX,0:NY+1), C(0:NX,0:NY+1)
      COMPLEX*16 G(0:NX,0:NY+1)
      COMPLEX*16 ICPACK(0:NX),OCPACK(0:NX)

      IF (RANK.NE.NPROCS-1) THEN
C If we aren't the highest process, then wait for data
        CALL MPI_RECV(OCPACK,NX+1,MPI_DOUBLE_COMPLEX,RANK+1,11
     &               ,MPI_COMM_WORLD,status,ierror)
        DO I=0,NX
          G(I,NY+1)=OCPACK(I)
        END DO
      ELSE
C Else, if we are the highest process, compute the solution directly at j=NY
        DO I=0,NX
          G(I,NY+1)=G(I,NY+1)/B(I,NY+1)
        END DO
      END IF

C All processes solve from NY-1..0
      DO J=NY,0,-1
      DO I=0,NX
        G(I,J)=(G(I,J)-C(I,J)*G(I,J+1))/B(I,J)
      END DO
      END DO

      IF (RANK.NE.0) THEN
        DO I=0,NX
          ICPACK(I)=G(I,2)
        END DO
        CALL MPI_SEND(ICPACK,NX+1,MPI_DOUBLE_COMPLEX,RANK-1,11
     &               ,MPI_COMM_WORLD,ierror)
      END IF

      RETURN
      END


!----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE INIT_MPI
C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
      INCLUDE 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'
      CHARACTER*35 MPI_IO_NUM 
! This subroutine initializes all mpi variables

      CALL MPI_INIT(IERROR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPROCS,IERROR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,RANK,IERROR)

      write(*,*) 'MPI Initialized, NPROCS,RANK: ',NPROCS,RANK

C Set a string to determine which input/output files to use
C When MPI is used, each process will read/write to files with the
C number of their rank (+1) appended to the end of the file.
C The string MPI_IO_NUM will be used to define the RANK+1 for each process
        IF (NPROCS.le.10) THEN
          MPI_IO_NUM=CHAR(MOD(RANK+1,10)+48)
        ELSE IF (NPROCS.le.100) THEN
          MPI_IO_NUM=CHAR(MOD(RANK+1,100)/10+48)
     &             //CHAR(MOD(RANK+1,10)+48)
        ELSE IF (NPROCS.le.1000) THEN
          MPI_IO_NUM=CHAR(MOD(RANK+1,1000)/100+48)
     &             //CHAR(MOD(RANK+1,100)/10+48)
     &             //CHAR(MOD(RANK+1,10)+48)
        ELSE IF (NPROCS.le.10000) THEN
          MPI_IO_NUM=CHAR(MOD(RANK+1,10000)/1000+48)
     &             //CHAR(MOD(RANK+1,1000)/100+48)
     &             //CHAR(MOD(RANK+1,100)/10+48)
     &             //CHAR(MOD(RANK+1,10)+48)
        ELSE
           WRITE(6,*) 'ERROR, NPROCS>10,000, Unsupported problem size'
        END IF

      RETURN
      END


!----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE INIT_CHAN_MPI
C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
C Initialize any constants here
      INCLUDE 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'

      INTEGER J, N

      if (rank .eq.0)then
      write(*,*) '*******IN INIT_CHAN_MPI *********'
      endif
 

      PI=4.D0*ATAN(1.D0)

! Set the upper and lower bounds for timestepping
      IF (RANK.eq.0) THEN
        write(*,*) 'U_BC_YMIN: ',U_BC_YMIN
        JEND=NY
        IF (U_BC_YMIN.EQ.0) THEN
          JSTART=2
        ELSE IF (U_BC_YMIN.EQ.1) THEN
          JSTART=1
        ELSE
          JSTART=2
        END IF
! Now, set the indexing for the scalar equations
        DO N=1,N_TH
          JEND_TH(N)=NY
          IF (TH_BC_YMIN(N).EQ.0) THEN
            JSTART_TH(N)=1
          ELSE IF (TH_BC_YMIN(N).EQ.1) THEN
            JSTART_TH(N)=0
          ELSE IF (TH_BC_YMIN(N).EQ.3) THEN
             IF (N ==1) THEN
              JSTART_TH(N)=1
             ELSE
              JSTART_TH(N)=0
             ENDIF
          ELSE
            JSTART_TH(N)=2
          END IF
        END DO
      ELSE IF (RANK.eq.NPROCS-1) THEN
        write(*,*) 'U_BC_YMAX: ',U_BC_YMAX
        JSTART=2
        IF (U_BC_YMAX.EQ.0) THEN
          JEND=NY-1
        ELSE IF (U_BC_YMAX.EQ.1) THEN
          JEND=NY
        ELSE
          JEND=NY-1
        END IF

! Set the upper and lower limits of timestepping of the scalar equations
        DO N=1,N_TH
        JSTART_TH(N)=2
        IF (TH_BC_YMAX(N).EQ.0) THEN
          JEND_TH(N)=NY
        ELSE IF (TH_BC_YMAX(N).EQ.1) THEN
          JEND_TH(N)=NY+1
        ELSE
          JEND_TH(N)=NY-1
        END IF
        END DO

      ELSE
! Here, we are on a middle process
        JSTART=2
        JEND=NY
        DO N=1,N_TH
          JSTART_TH(N)=2
          JEND_TH(N)=NY
        ENDDO
      END IF

      RETURN
      END

      SUBROUTINE APPLY_BC_TH_MPI(MATL,MATD,MATU,VEC,K,N)
! This subroutine applies the boundary conditions to the
! scalar fields prior to the implicit solve
      INCLUDE 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi' 
      INTEGER k, N
! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
        IF (RANK.eq.0) THEN
! If we have the lowest plane, apply the boundary conditions
          CALL APPLY_BC_TH_LOWER(MATL,MATD,MATU,VEC,K,N)
        END IF
        IF (RANK.eq.NPROCS-1) THEN
! If we have the upper plane, apply the boundary conditions
          CALL APPLY_BC_TH_UPPER(MATL,MATD,MATU,VEC,N)
        END IF
      RETURN
      END

      SUBROUTINE APPLY_BC_TH_MPI_C(MATL_C,MATD_C,MATU_C,VEC_C,N)
! This subroutine applies the boundary conditions to the
! scalar fields prior to the implicit solve
      INCLUDE 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'
      INTEGER N
! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
        IF (RANK.eq.0) THEN
! If we have the lowest plane, apply the boundary conditions
          CALL APPLY_BC_TH_LOWER_C(MATL_C,MATD_C,MATU_C,VEC_C,N)
        END IF
        IF (RANK.eq.NPROCS-1) THEN
! If we have the upper plane, apply the boundary conditions
          CALL APPLY_BC_TH_UPPER_C(MATL_C,MATD_C,MATU_C,VEC_C,N)
        END IF
      RETURN
      END


      SUBROUTINE APPLY_BC_U2_MPI(MATL,MATD,MATU,VEC,K)
! This subroutine applies the boundary conditions to the
! velocity field prior to the implicit solve
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'
      integer k
! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
        IF (RANK.eq.0) THEN
! If we have the lowest plane, apply the boundary conditions
          CALL APPLY_BC_2_LOWER(MATL,MATD,MATU,VEC,k)
        END IF
        IF (RANK.eq.NPROCS-1) THEN
! If we have the highest plane, apply the boundary conditions 
          CALL APPLY_BC_2_UPPER(MATL,MATD,MATU,VEC,k)
        END IF
      RETURN
      END

      SUBROUTINE APPLY_BC_U2_MPI_C(MATL_C,MATD_C,MATU_C,VEC_C)
! This subroutine applies the boundary conditions to the
! velocity field prior to the implicit solve
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'
! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
        IF (RANK.eq.0) THEN
! If we have the lowest plane, apply the boundary conditions
          CALL APPLY_BC_2_LOWER_C(MATL_C,MATD_C,MATU_C,VEC_C)
        END IF
        IF (RANK.eq.NPROCS-1) THEN
! If we have the highest plane, apply the boundary conditions 
          CALL APPLY_BC_2_UPPER_C(MATL_C,MATD_C,MATU_C,VEC_C)
        END IF
      RETURN
      END


      SUBROUTINE APPLY_BC_U1_MPI(MATL,MATD,MATU,VEC,k)
! This subroutine applies the boundary conditions to the
! velocity field prior to the implicit solve
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'
      integer I,J,k
! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
        IF (RANK.eq.0) THEN
! If we have the lowest plane, apply the boundary conditions
          CALL APPLY_BC_1_LOWER(MATL,MATD,MATU,VEC,k)        
        END IF
        IF (RANK.eq.NPROCS-1) THEN
! If we have the highest plane, apply the boundary conditions 
          CALL APPLY_BC_1_UPPER(MATL,MATD,MATU,VEC,k)
        END IF
      RETURN
      END

      SUBROUTINE APPLY_BC_U1_MPI_C(MATL_C,MATD_C,MATU_C,VEC_C)
! This subroutine applies the boundary conditions to the
! velocity field prior to the implicit solve
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'
! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
        IF (RANK.eq.0) THEN
! If we have the lowest plane, apply the boundary conditions
          CALL APPLY_BC_1_LOWER_C(MATL_C,MATD_C,MATU_C,VEC_C)
        END IF
        IF (RANK.eq.NPROCS-1) THEN
! If we have the highest plane, apply the boundary conditions 
          CALL APPLY_BC_1_UPPER_C(MATL_C,MATD_C,MATU_C,VEC_C)
        END IF
      RETURN
      END


      SUBROUTINE APPLY_BC_U3_MPI(MATL,MATD,MATU,VEC,k)
! This subroutine applies the boundary conditions to the
! velocity field prior to the implicit solve
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi' 
      integer k
! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
        IF (RANK.eq.0) THEN
! If we have the lowest plane, apply the boundary conditions
          CALL APPLY_BC_3_LOWER(MATL,MATD,MATU,VEC,k)
        END IF
        IF (RANK.eq.NPROCS-1) THEN
! If we have the highest plane, apply the boundary conditions 
          CALL APPLY_BC_3_UPPER(MATL,MATD,MATU,VEC,k)
        END IF
      RETURN
      END

      SUBROUTINE APPLY_BC_U3_MPI_C(MATL_C,MATD_C,MATU_C,VEC_C)
! This subroutine applies the boundary conditions to the
! velocity field prior to the implicit solve
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'
! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
        IF (RANK.eq.0) THEN
! If we have the lowest plane, apply the boundary conditions
          CALL APPLY_BC_3_LOWER_C(MATL_C,MATD_C,MATU_C,VEC_C)
        END IF
        IF (RANK.eq.NPROCS-1) THEN
! If we have the highest plane, apply the boundary conditions 
          CALL APPLY_BC_3_UPPER_C(MATL_C,MATD_C,MATU_C,VEC_C)
        END IF
      RETURN
      END


      SUBROUTINE APPLY_BC_REM_DIV_MPI(MATL_C,MATD_C,MATU_C,VEC_C,K)
! This subroutine applies the boundary conditions for the Poisson Eq.
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'
      INTEGER I,J,K
! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
C Apply the boundary conditions
        IF (RANK.EQ.0) THEN
        DO I=0,NKX
C Use homogeneous dirichlet BCS for kx=kz=0 component at bottom wall
          IF ((K.EQ.0).AND.(I.EQ.0)) THEN
C Otherwise the matrix will be singular
C Use homogeneous dirichlet BCS for kx=kz=0 component at bottom wall
            MATL_C(I,1)=0.
            MATD_C(I,1)=1.
            MATU_C(I,1)=0.
            VEC_C(I,1)=(0.,0.)
          ELSE
C Use Dirichlet boundary conditions, dp/dz=0 at walls
            MATL_C(I,1)=0.
            MATD_C(I,1)=1.
            MATU_C(I,1)=-1.
            VEC_C(I,1)=(0.,0.)
          END IF
        END DO
        END IF
C Apply the boundary conditions
        IF (RANK.EQ.NPROCS-1) THEN
          DO I=0,NKX
            MATL_C(I,NY)=1.
            MATD_C(I,NY)=-1.
            MATU_C(I,NY)=0.
            VEC_C(I,NY)=(0.,0.)
          END DO
        END IF
        RETURN
        END

      SUBROUTINE APPLY_BC_POISSON_MPI(MATL_C,MATD_C,MATU_C,VEC_C,K)
! This subroutine applies the boundary conditions for the Poisson Eq.
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'
      INTEGER I,K
! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
C Use dirichlet boundary condition at the lower wall to
C prevent the tridiagonal matrix from becomming singular for i,k=0
        IF (RANK.eq.0) THEN
        DO I=0,NKX
          IF ((I.EQ.0).AND.(K.EQ.0)) THEN
            MATD_C(I,1)=1.
            MATU_C(I,1)=0.
            VEC_C(I,1)=(0.,0.)
          ELSE
! Here, apply Neumann boundary conditions (dp/dz=0) at the walls
            MATD_C(I,1)=1.
            MATU_C(I,1)=-1.
            VEC_C(I,1)=(0.,0.)
          END IF
        END DO
        END IF
C Use dirichlet boundary condition at the lower wall to
C prevent the tridiagonal matrix from becomming singular for i,k=0
        IF (RANK.eq.NPROCS-1) THEN
        DO I=0,NKX
            MATD_C(I,NY)=-1.
            MATL_C(I,NY)=1.
            VEC_C(I,NY)=(0.,0.)
        END DO
        END IF

      RETURN
      END

      SUBROUTINE APPLY_WALL_MODEL_MPI
! This subroutine applies the boundary conditions for the Poisson Eq.
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'

! Apply Boundary conditions to velocity field
      IF (RANK.EQ.0) THEN
        if ((U_BC_YMIN.EQ.3).or.(W_BC_YMIN.EQ.3)) then
         CALL WALL_MODEL_LOWER
        endif
      END IF

      IF (RANK.EQ.NPROCS-1) THEN
        if ((U_BC_YMAX.EQ.3).or.(W_BC_YMAX.EQ.3)) then
         CALL WALL_MODEL_UPPER
        endif
      END IF

      RETURN
      END

      SUBROUTINE APPLY_BC_VEL_MPI
! This subroutine applies the boundary conditions for the Poisson Eq.
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'

! Apply Boundary conditions to velocity field
      IF (RANK.EQ.0) THEN
        J1i=2
      IF ((U_BC_YMIN.EQ.3).or.(W_BC_YMIN.EQ.3)) THEN
        CALL APPLY_BC_NWM_LOWER
        J1 = 2
      ELSE  
c        J1 = JSTART
c        J1 = 1
c        J1i=2
        CALL APPLY_BC_VEL_LOWER
      END IF
      END IF

      IF (RANK.EQ.NPROCS-1) THEN
        J2e=NY
      IF ((U_BC_YMAX.EQ.3).or.(W_BC_YMAX.EQ.3)) then
        CALL APPLY_BC_NWM_UPPER
        J2 = NY-1
      ELSE
        CALL APPLY_BC_VEL_UPPER
c        J2 = JEND
c        J2 = NY 
      END IF
      END IF 
 
      RETURN
      END
 


      SUBROUTINE APPLY_BC_ANWM_MPI
! This subroutine applies the boundary conditions for the Poisson Eq.
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'

! Apply Boundary conditions to velocity field
      IF (RANK.EQ.0) THEN
        CALL APPLY_BC_ANWM_LOWER
      END IF

      IF (RANK.EQ.NPROCS-1) THEN
        CALL APPLY_BC_ANWM_UPPER
      END IF

      RETURN
      END


      SUBROUTINE APPLY_NUT_NWM_MPI
      INCLUDE 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'
      integer i,k,N

       IF (RANK.EQ.0) THEN
            DO I=0,NXM
            DO K=0,NZM
              NU_T(I,K,2)=0.d0
              DO N=1,N_TH
                KAPPA_T(I,K,2,N)=0.d0
              END DO
            END DO
            END DO
       END IF

       IF (RANK.EQ.NPROCS-1) THEN
            DO I=0,NXM
            DO K=0,NZM
              NU_T(I,K,NY)=0.d0
              DO N=1,N_TH
                KAPPA_T(I,K,NY,N)=0.d0
              END DO
            END DO
            END DO
       END IF

       RETURN
       END 


      SUBROUTINE APPLY_NUT_NWM_MPI_UPPER
      INCLUDE 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'
      integer i,k,N


       IF (RANK.EQ.NPROCS-1) THEN
            DO I=0,NXM
            DO K=0,NZM
              NU_T(I,K,NY)=0.d0
              DO N=1,N_TH
                KAPPA_T(I,K,NY,N)=0.d0
              END DO
            END DO
            END DO
       END IF

       RETURN
       END



      SUBROUTINE APPLY_NUT_NWM_MPI_LOWER
      INCLUDE 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'
      integer i,k,N

       IF (RANK.EQ.0) THEN
            DO I=0,NXM
            DO K=0,NZM
              NU_T(I,K,2)=0.d0
              DO N=1,N_TH
                KAPPA_T(I,K,2,N)=0.d0
              END DO
            END DO
            END DO
       END IF


       RETURN
       END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE THOMAS_FORWARD_REAL_MPI_BACK(A,B,C,G,NY,NX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This subroutine performs the backward sweep of the Thomas algorithm
C Thomas algorithm solves Ax=b for tridiagonal A
C The RHS vector and solution are real
C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
C Returns solution in x
C The indexing should be done by ROW, ie.
C [ b1  c1   0   0   0 ...
C [ a2  b2  c2   0   0 ...
C [  0  a3  b3   c3  0 ...

      INCLUDE 'mpif.h'
      INCLUDE 'header_mpi'

      INTEGER I,J,NY
      REAL*8 A(0:NX-1,0:NY+1), B(0:NX-1,0:NY+1)
     &  , C(0:NX-1,0:NY+1), G(0:NX-1,0:NY+1)

      REAL*8 OCPACK(3),ICPACK(3)

      DO I=0,NX-1

      IF (RANK.NE.0) THEN
C If we aren't the lowest process, then wait for data
        CALL MPI_RECV(OCPACK,3,MPI_DOUBLE_PRECISION,RANK-1,12
     &               ,MPI_COMM_WORLD,status,ierror)
C Unpack the data
        B(I,1)=OCPACK(1)
        C(I,1)=OCPACK(2)
        G(I,1)=OCPACK(3)
      ELSE
C If we are the lowest process, start at J=1
        J=1
        A(I,J)=-A(I,J)/B(I,J-1)
        B(I,J)=B(I,J)+A(I,J)*C(I,J-1)
        G(I,J)=G(I,J)+A(I,J)*G(I,J-1)
      END IF

      DO J=2,NY
        A(I,J)=-A(I,J)/B(I,J-1)
        B(I,J)=B(I,J)+A(I,J)*C(I,J-1)
        G(I,J)=G(I,J)+A(I,J)*G(I,J-1)
      END DO

      IF (RANK.NE.NPROCS-1) THEN
        ICPACK(1)=B(I,NY)
        ICPACK(2)=C(I,NY)
        ICPACK(3)=G(I,NY)
        CALL MPI_SEND(ICPACK,3,MPI_DOUBLE_PRECISION,RANK+1,12
     &               ,MPI_COMM_WORLD,ierror)
      END IF

      END DO

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE THOMAS_BACKWARD_REAL_MPI_BACK(A,B,C,G,NY,NX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This subroutine performs the backward sweep of the Thomas algorithm
C Thomas algorithm solves Ax=b for tridiagonal A
C The RHS vector and solution are real
C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
C Returns solution in x
C The indexing should be done by ROW, ie.
C [ b1  c1   0   0   0 ...
C [ a2  b2  c2   0   0 ...
C [  0  a3  b3   c3  0 ...

      INCLUDE 'mpif.h'
      INCLUDE 'header_mpi'

      INTEGER I,J,NY
      REAL*8 A(0:NX-1,0:NY+1), B(0:NX-1,0:NY+1), C(0:NX-1,0:NY+1)
     &         , G(0:NX-1,0:NY+1)

      DO I=0,NX-1

      IF (RANK.NE.NPROCS-1) THEN
C If we aren't the highest process, then wait for data
        CALL MPI_RECV(G(I,NY+1),1,MPI_DOUBLE_PRECISION,RANK+1,13
     &               ,MPI_COMM_WORLD,status,ierror)
      ELSE
C Else, if we are the highest process, then compute the solution at j=NY
        G(I,NY)=G(I,NY)/B(I,NY)
      END IF

      DO J=NY,1,-1
        G(I,J)=(G(I,J)-C(I,J)*G(I,J+1))/B(I,J)
      END DO

      IF (RANK.NE.0) THEN
        CALL MPI_SEND(G(I,1),1,MPI_DOUBLE_PRECISION,RANK-1,13
     &               ,MPI_COMM_WORLD,ierror)
      END IF

      END DO

      RETURN
      END

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE MPI_INST_PRO(OPACK,IPACK)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      include 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'
      integer :: send_rank!,en,st
      REAL(8),DIMENSION(0:NXM,1:NY*NPROCS) :: IPACK
      REAL(8),DIMENSION(0:NXM,1:NY) :: OPACK
!      REAL(8),DIMENSION(0:NXM,1:NY,1:NPROCS) ::temp

      send_rank=0
! READ IN LOCAL DT BASED ON CFL
!      OPACK = DT
      
! SEND ALL LOCAL DT'S TO EVERY PROCESS
      CALL MPI_GATHER(OPACK,NY*(NXM+1),MPI_DOUBLE_PRECISION,IPACK,NY
     + *(NXM+1), MPI_DOUBLE_PRECISION,send_rank,MPI_COMM_WORLD,IERROR)

!     T = MINVAL(IPACK)
!       IF ( RANK .EQ. 0) THEN
!        temp(:,:) = 0.0d0
!        st=1
!        DO I = 1,NPROCS
!          en=st+SIZE_1-1
!          DO J =1,SIZE_2
!           temp(st:en,J) =  STAT_TMP(:,J,I)
!          ENDDO
!          st=en+1
!        ENDDO
!      ENDIF

      RETURN
      END

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE MPI_PROF_ALL_STAT
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      include 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'



      call MPI_INST_PROF(dble(CR1(0,0,1:NY)),UME_TOT)
      call MPI_INST_PROF(dble(CR2(0,0,1:NY)),VME_TOT)
      call MPI_INST_PROF(dble(CR3(0,0,1:NY)),WME_TOT)
      call MPI_INST_PROF(dble(CRTH(0,0,1:NY,1)),TH_TOT)
      call MPI_INST_PROF(dble(CRTH(0,0,1:NY,2)),SAL_TOT)
      call MPI_INST_PROF(urms(1:NY),urms_tot)
      call MPI_INST_PROF(vrms(1:NY),vrms_tot)
      call MPI_INST_PROF(wrms(1:NY),wrms_tot)
      call MPI_INST_PROF(shear(1:NY),shear_t)
      call MPI_INST_PROF(dthdyt(1:NY),dthdy_t)
      call MPI_INST_PROF(dthdys(1:NY),dsdy_t)
      call MPI_INST_PROF(dudy(1:NY),dudy_t)
      call MPI_INST_PROF(dwdy(1:NY),dwdy_t)
      call MPI_INST_PROF(thrms(1:NY,1),thrms_t)
      call MPI_INST_PROF(thrms(1:NY,2),srms_t)
      call MPI_INST_PROF(thus(1:NY),thus_t)
      call MPI_INST_PROF(thut(1:NY),thut_t)
      call MPI_INST_PROF(thvs(1:NY),thvs_t)
      call MPI_INST_PROF(thvt(1:NY),thvt_t)
      call MPI_INST_PROF(uv(1:NY),uv_t)
      call MPI_INST_PROF(wv(1:NY),wv_t)
      call MPI_INST_PROF(uw(1:NY),uw_t)
      call MPI_INST_PROF(tke_mean(1:NY),tke_t)
      call MPI_INST_PROF(tke_7(1:NY),diss_t)
      call MPI_INST_PROF(tke_3(1:NY),prod_t)
      call MPI_INST_PROF(bu_pro(1:NY),bpro_t)

      RETURN
      END
