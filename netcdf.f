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

      dim_2d_t(1) = y_dimid
      dim_2d_t(2) = z_dimid
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
      retval=nf_def_var(ncid_stats,"Temp",NF_DOUBLE,2,dim_1d_t
     &        ,th_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_stats,"Sal",NF_DOUBLE,2,dim_1d_t
     &        ,sal_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_stats,"Temp_rms",NF_DOUBLE,2,dim_1d_t
     &        ,thr_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_stats,"Sal_rms",NF_DOUBLE,2,dim_1d_t
     &        ,salr_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_stats,"Temp_grad",NF_DOUBLE,2,dim_1d_t
     &        ,dthdy_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_stats,"Sal_grad",NF_DOUBLE,2,dim_1d_t
     &        ,dsdy_varid)
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
      retval=nf_def_var(ncid_stats,"dudy",NF_DOUBLE,2,dim_1d_t
     &        ,dudy_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_stats,"dwdy",NF_DOUBLE,2,dim_1d_t
     &        ,dwdy_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
       retval=nf_def_var(ncid_stats,"ut_prime",NF_DOUBLE,2,dim_1d_t
     &        ,thut_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_stats,"us_prime",NF_DOUBLE,2,dim_1d_t
     &        ,thus_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_stats,"wt_prime",NF_DOUBLE,2,dim_1d_t
     &        ,thvt_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_stats,"ws_prime",NF_DOUBLE,2,dim_1d_t
     &        ,thvs_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_stats,"tke",NF_DOUBLE,2,dim_1d_t
     &        ,tke_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_stats,"buo_prod",NF_DOUBLE,2,dim_1d_t
     &        ,bp_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_stats,"shear_prod",NF_DOUBLE,2,dim_1d_t
     &        ,sp_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_stats,"dissip",NF_DOUBLE,2,dim_1d_t
     &        ,di_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
       retval=nf_def_var(ncid_stats,"U2",NF_DOUBLE,3,dim_2d_t
     &        ,u2_varid)
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
!      nc_count(1)=NY_TOT
!      nc_count(2)=1
!      nc_start(1)=1
       nc_start=(/1, 1/)
       nc_count=(/NY_TOT, 1/)
!counts for the 2d variables
       nc_start2=(/1, 1, 1/)
       nc_count2=(/NY_TOT, NZ, 1/)
       nc_index=1
      endif
      RETURN
      END


      subroutine NETCDF_WRITE_STATS_CHAN
      include 'header'
      include 'netcdf.inc'
       integer retval,j,dum
       !dum=61960+140
       write(6,*) 'writing netcdf'
!      nc_index=nc_index+1
       dum=dum
!      write(6,*)TIME_STEP,dum,TIME_STEP-dum
!      nc_start(2)=(nc_index+nc_start(1))!TIME_STEP-dum)/2
!      nc_start2(3)=(nc_index+nc_start2(1))
!       nc_start(2)=TIME_STEP-dum+1
       write(6,*) nc_start(1),nc_start(2)
C Write the data to the netCDF file
       retval=nf_sync(ncid_stats)
C      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid_stats,ubar_varid
     &             ,nc_start,nc_count,UME_TOT(1:NY_TOT))
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid_stats,vbar_varid
     &             ,nc_start,nc_count,VME_TOT(1:NY_TOT))
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid_stats,wbar_varid
     &             ,nc_start,nc_count,WME_TOT(1:NY_TOT))
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid_stats,urms_varid
     &             ,nc_start,nc_count,urms_tot(1:NY_TOT))
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
      retval = nf_put_vara_double(ncid_stats,wv_varid
     &             ,nc_start,nc_count,wv_t(1:NY_TOT))
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid_stats,uw_varid
     &             ,nc_start,nc_count,uw_t(1:NY_TOT))
      retval = nf_put_vara_double(ncid_stats,th_varid
     &             ,nc_start,nc_count,TH_TOT(1:NY_TOT))
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid_stats,sal_varid
     &             ,nc_start,nc_count,SAL_TOT(1:NY_TOT))
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid_stats,thr_varid
     &             ,nc_start,nc_count,thrms_t(1:NY_TOT))
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid_stats,salr_varid
     &             ,nc_start,nc_count,srms_t(1:NY_TOT))
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid_stats,dthdy_varid
     &             ,nc_start,nc_count,dthdy_t(1:NY_TOT))
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid_stats,dsdy_varid
     &             ,nc_start,nc_count,dsdy_t(1:NY_TOT))
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid_stats,dudy_varid
     &             ,nc_start,nc_count,dudy_t(1:NY_TOT))
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid_stats,dwdy_varid
     &             ,nc_start,nc_count,dwdy_t(1:NY_TOT))
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid_stats,thut_varid
     &             ,nc_start,nc_count,thut_t(1:NY_TOT))
      retval = nf_put_vara_double(ncid_stats,thus_varid
     &             ,nc_start,nc_count,thus_t(1:NY_TOT))
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid_stats,thvt_varid
     &             ,nc_start,nc_count,thvt_t(1:NY_TOT))
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid_stats,thvs_varid
     &             ,nc_start,nc_count,thvs_t(1:NY_TOT))
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid_stats,tke_varid
     &             ,nc_start,nc_count,tke_t(1:NY_TOT))
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid_stats,bp_varid
     &             ,nc_start,nc_count,bpro_t(1:NY_TOT))
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid_stats,sp_varid
     &             ,nc_start,nc_count,prod_t(1:NY_TOT))
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid_stats,di_varid
     &             ,nc_start,nc_count,diss_t(1:NY_TOT))
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid_stats,u2_varid
     &             ,nc_start,nc_count,U2_ty(1:NY_TOT,0:NZ-1))
      if (retval .ne. nf_noerr) call handle_err(retval)
      
C Write the current time to the netCDF file
      retval = nf_put_vara_double(ncid_stats,t_varid,
     &             nc_index, 1, TIME)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_sync(ncid_stats)
      if (retval .ne. nf_noerr) call handle_err(retval)
      write(6,*)'syncrhonized data'
      nc_index=nc_index+1
      nc_start(2)=nc_start(2)+1
      nc_start2(3)=nc_start2(3)+1
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
      integer retval,len
      write(6,*) 'I am reading netcdf now' ,int(TIME)
      retval =nf_open('mean_prof_stat.nc',NF_WRITE,ncid_stats)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval =nf_inq_dimid(ncid_stats,"t",t_dimid)      
      if (retval .ne. nf_noerr) call handle_err(retval)
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
      retval =nf_inq_varid(ncid_stats,"Temp",th_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval =nf_inq_varid(ncid_stats,"Sal",sal_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval =nf_inq_varid(ncid_stats,"Temp_rms",thr_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval =nf_inq_varid(ncid_stats,"Sal_rms",salr_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval =nf_inq_varid(ncid_stats,"Temp_grad",dthdy_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval =nf_inq_varid(ncid_stats,"Sal_grad",dsdy_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval =nf_inq_varid(ncid_stats,"dudy",dudy_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval =nf_inq_varid(ncid_stats,"dwdy",dwdy_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval =nf_inq_varid(ncid_stats,"ut_prime",thut_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval =nf_inq_varid(ncid_stats,"us_prime",thus_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval =nf_inq_varid(ncid_stats,"wt_prime",thvt_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval =nf_inq_varid(ncid_stats,"ws_prime",thvs_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval =nf_inq_varid(ncid_stats,"tke",tke_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval =nf_inq_varid(ncid_stats,"buo_prod",bp_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval =nf_inq_varid(ncid_stats,"shear_prod",sp_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval =nf_inq_varid(ncid_stats,"dissip",di_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval =nf_inq_varid(ncid_stats,"U2",u2_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval =nf_inq_dimlen(ncid_stats,t_dimid,len)
!      nc_count(1)=NY_TOT
!      nc_count(2)=1
!      nc_start(1)=1
       nc_start(2)=len

!      nc_count2(1)=NY_TOT
!      nc_count2(2)=NZ
!      nc_count2(3)=1
!      nc_start2(1)=1
!      nc_start2(2)=1
       nc_start2(3)=len

!      nc_index=len-1
      write(6,*) 'i am here so far',len
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
