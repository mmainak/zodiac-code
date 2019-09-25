      subroutine courant
! This subroutine sets the timestep based on the specified CFL number
! The subroutine should be called with the velocity in physical space

      include 'header'

      real*8 vel 

      dt=0.d0

      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            dt_x=cfl*dx(i)/U1(i,k,j)
            dt_y=cfl*dy(j)/U2(i,k,j)         
            dt_z=cfl*dz(k)/U3(i,k,j)
            dt=min(dt,dt_x,dt_y,dt_z)
          end do
        end do
      end do 

      write(*,*) 'CFL #, DT: ',cfl,dt

      if (dt.le.0.d0) then
        pause 'Error: dt<=0 in courant: CFL number is too restrictive'
      else
        DELTA_T=dt
        H_BAR(1)=DELTA_T*(8.0/15.0)
        H_BAR(2)=DELTA_T*(2.0/15.0)
        H_BAR(3)=DELTA_T*(5.0/15.0)
      end if 

      return
      end

  
