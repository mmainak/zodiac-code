C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE NETCDF_OPEN_STATS_CHAN
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INCLUDE 'netcdf.inc'

      integer dim_2d(2)
      integer dim_2d_t(3), dim_1d_t(2)
      integer retval
      logical there

      integer i,j
!Inquire if the netcdf file exists from previous run. Otherwise create a
!new netcdf file
      inquire(file='mean_prof_stat.nc',exist=there)
      if (there) then 
       write(6,*) 'checking netcdf'
       call NETCDF_READ_STATS
      else
! Create a new netCDF file    
 
      retval = nf_create('mean_prof_stat.nc', NF_NOCLOBBER
     &                   , ncid_stats)
      if (retval .ne. nf_noerr) call handle_err(retval)


!Define the dimensions. NetCDF will hand back an ID for each.
      retval = nf_def_dim(ncid_stats, "x", NX, x_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid_stats, "y", NY_TOT, y_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid_stats, "z", NZ, z_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid_stats, "t", NF_UNLIMITED, t_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)

! Define the coordinate variables in the netCDF file
      retval=nf_def_var(ncid_stats,"x",NF_DOUBLE,1,x_dimid,x_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_stats,"y",NF_DOUBLE,1,y_dimid,y_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_stats,"z",NF_DOUBLE,1,z_dimid,z_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_stats,"t",NF_DOUBLE,1,t_dimid,t_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

!     The dimids array is used to pass the IDs of the dimensions of
!     the variables. Note that in fortran arrays are stored in
!     column-major format.

      dim_2d(1) = x_dimid
      dim_2d(2) = y_dimid

      dim_2d_t(1) = x_dimid
      dim_2d_t(2) = y_dimid
      dim_2d_t(3) = t_dimid

      dim_1d_t(1) = y_dimid
      dim_1d_t(2) = t_dimid

! Define the variable data in our netCDF file
      retval=nf_def_var(ncid_stats,"U_BAR",NF_DOUBLE,2,dim_1d_t
     &        ,ubar_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_stats,"V_BAR",NF_DOUBLE,2,dim_1d_t
     &        ,vbar_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_stats,"W_BAR",NF_DOUBLE,2,dim_1d_t
     &        ,wbar_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_stats,"U_RMS",NF_DOUBLE,2,dim_1d_t
     &        ,urms_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_stats,"V_RMS",NF_DOUBLE,2,dim_1d_t
     &        ,vrms_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_stats,"W_RMS",NF_DOUBLE,2,dim_1d_t
     &        ,wrms_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_stats,"uv",NF_DOUBLE,2,dim_1d_t
     &        ,uv_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_stats,"wv",NF_DOUBLE,2,dim_1d_t
     &        ,wv_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_stats,"uw",NF_DOUBLE,2,dim_1d_t
     &        ,uw_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

!End define mode. This tells netCDF we are done defining metadata.
      retval = nf_enddef(ncid_stats)
      if (retval .ne. nf_noerr) call handle_err(retval)

!Now put the grid data 
    
      retval = nf_put_var_double(ncid_stats,x_varid,GX(0:NX-1))
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_double(ncid_stats,z_varid,GZ(0:NZ-1))
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_double(ncid_stats,y_varid,GYF_TOT(1:NY_TOT))
      if (retval .ne. nf_noerr) call handle_err(retval)

! Define the start and count arrays
      nc_count(1)=NY_TOT
      nc_count(2)=1
      nc_start(1)=1
      nc_index=0
      endif
      RETURN
      END


      subroutine NETCDF_WRITE_STATS_CHAN
      include 'header'
      include 'netcdf.inc'
       integer retval,j,dum
       dum=61960
       write(6,*) 'writing netcdf'
       nc_index=nc_index+1
       dum=dum
       write(6,*)TIME_STEP-dum
      nc_start(2)=nc_index+nc_start(1)+TIME_STEP-dum
       write(6,*) nc_start(1),nc_start(2)
C Write the data to the netCDF file
C       retval=nf_sync(ncid_stats)
C      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid_stats,ubar_varid
     &             ,nc_start,nc_count,UME_TOT(1:NY_TOT))
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid_stats,ubar_varid
     &             ,nc_start,nc_count,UME_TOT(1:NY_TOT))
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid_stats,vbar_varid
     &             ,nc_start,nc_count,VME_TOT(1:NY_TOT))
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid_stats,wbar_varid
     &             ,nc_start,nc_count,WME_TOT(1:NY_TOT))
      if (retval .ne. nf_noerr) call handle_err(retval)
       
      retval = nf_put_vara_double(ncid_stats,vrms_varid
     &             ,nc_start,nc_count,vrms_tot(1:NY_TOT))
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid_stats,wrms_varid
     &             ,nc_start,nc_count,wrms_tot(1:NY_TOT))
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid_stats,uv_varid
     &             ,nc_start,nc_count,uv_t(1:NY_TOT))
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid_stats,uw_varid
     &             ,nc_start,nc_count,uw_t(1:NY_TOT))
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid_stats,wv_varid
     &             ,nc_start,nc_count,wv_t(1:NY_TOT))
      
C Write the current time to the netCDF file
      retval = nf_put_vara_double(ncid_stats,t_varid,
     &             nc_index, 1, TIME)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_sync(ncid_stats)
      if (retval .ne. nf_noerr) call handle_err(retval)
      write(6,*)'syncrhonized data'
!      retval = nf_close(ncid_stats)
!      if (retval .ne. nf_noerr) call handle_err(retval)
      return
      end


      subroutine NETCDF_CLOSE_STATS_CHAN
      include 'header'
      include 'netcdf.inc'
      integer retval

C     Close the file. This frees up any internal netCDF resources
C     associated with the file, and flushes any buffers.
      retval = nf_close(ncid_stats)
      if (retval .ne. nf_noerr) call handle_err(retval)
    
      return
      end

      subroutine NETCDF_READ_STATS
      include 'header'
      include 'netcdf.inc'
      integer retval
      write(6,*) 'I am reading netcdf now' ,int(TIME)
      retval =nf_open('mean_prof_stat.nc',NF_WRITE,ncid_stats)
      if (retval .ne. nf_noerr) call handle_err(retval)
!      retval=nf_sync(ncid_stats)
!      if (retval .ne. nf_noerr) call handle_err(retval)
!      retval =nf_inq_dimid(ncid_stats,"x",x_dimid)
!      if (retval .ne. nf_noerr) call handle_err(retval)
!      retval =nf_inq_dimid(ncid_stats,"y",y_dimid)
!      if (retval .ne. nf_noerr) call handle_err(retval)
!      retval =nf_inq_dimid(ncid_stats,"z",z_dimid)
!      if (retval .ne. nf_noerr) call handle_err(retval)
!      retval =nf_inq_varid(ncid_stats,"t",t_dimid)      
!      if (retval .ne. nf_noerr) call handle_err(retval)
      retval =nf_inq_varid(ncid_stats,"x",x_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval =nf_inq_varid(ncid_stats,"y",y_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval =nf_inq_varid(ncid_stats,"z",z_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval =nf_inq_varid(ncid_stats,"t",t_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval =nf_inq_varid(ncid_stats,"U_BAR",ubar_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval =nf_inq_varid(ncid_stats,"V_BAR",vbar_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval =nf_inq_varid(ncid_stats,"W_BAR",wbar_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval =nf_inq_varid(ncid_stats,"U_RMS",urms_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval =nf_inq_varid(ncid_stats,"V_RMS",vrms_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval =nf_inq_varid(ncid_stats,"W_RMS",wrms_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval =nf_inq_varid(ncid_stats,"uv",uv_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval =nf_inq_varid(ncid_stats,"wv",wv_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval =nf_inq_varid(ncid_stats,"uw",uw_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
!     write(6,*) 'i am here so far' ,nc_count(1),nc_count(2)
      nc_count(1)=NY_TOT
      nc_count(2)=1
      nc_start(1)=1
      nc_index=0
      write(6,*) 'i am here so far'
      retval=nf_sync(ncid_stats)
      if (retval .ne. nf_noerr) call handle_err(retval)
      return

      end


      subroutine handle_err(errcode)
       include 'header'
       include 'netcdf.inc'
       integer, intent ( in) :: errcode

        if(errcode /= nf_noerr) then
           write(6,*) 'somethings wrong'
           print *, trim(nf_strerror(errcode))
           stop "Stopped"
        end if

      return
      end
