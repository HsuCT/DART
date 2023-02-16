program convert_gdc_cdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use          types_mod, only : r8
use   time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, set_time,&
                               increment_time, get_time, set_date, operator(-),  &
                               print_date
use      utilities_mod, only : initialize_utilities, find_namelist_in_file,    &
                               check_namelist_read, nmlfileunit, do_nml_file,  &
                               get_next_filename, error_handler, E_ERR, E_MSG, &
                               nc_check, find_textfile_dims, do_nml_term,      &
                               finalize_utilities
use       location_mod, only : VERTISHEIGHT, set_location
use   obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq,         &
                               static_init_obs_sequence, init_obs, destroy_obs,   &
                               write_obs_seq, init_obs_sequence,                  &
                               destroy_obs_sequence, set_obs_values, set_obs_def, &
                               set_copy_meta_data, set_qc, set_qc_meta_data
use   obs_def_mod,      only : obs_def_type, set_obs_def_time, set_obs_def_type_of_obs, &
                               set_obs_def_error_variance, set_obs_def_location,        &
                               set_obs_def_key
use       obs_kind_mod, only : GDC_VELOCITY_U, GDC_VELOCITY_V, GDC_TEMPERATURE, &
                               GDC_DENSITY, GDC_TEMPERATURE_ION, GDC_DENSITY_ION_OP, &
                               GDC_DENSITY_NEUTRAL_O, GDC_DENSITY_NEUTRAL_O2, GDC_DENSITY_NEUTRAL_N2
use  obs_utilities_mod, only : add_obs_to_seq

use           netcdf
implicit none
! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = &
   "$URL: https://svn-dares-dart.cgd.ucar.edu/DART/branches/recam/observations/obs_converters/mighti/convert_gdc_cdf.f90 $"
character(len=*), parameter :: revision = "$Revision: 12384 $"
character(len=*), parameter :: revdate  = "$Date: 2019-12-04 $"


integer, parameter ::   num_copies = 1,   &   ! number of copies in sequence
                        num_qc     = 1        ! number of QC entries

character (len=512) :: msgstring, msgstring2, next_infile
character (len=19)  :: datestr
character (len=NF90_MAX_NAME)  :: name
integer :: ncid, varid, nlevels, k,m, nfiles, num_new_obs, oday, osec, aday, asec, &
           iyear, imonth, iday, ihour, imin, isec, zloc, obs_num, &
           io, iunit, nobs, filenum, dummy, numrejected, st_sec, st_day,   &
           timestamps,vec_comp,horizontal
logical :: file_exist, first_obs, from_list = .false.
real(r8) :: oerr, lato, lono, hghto, obs_window,     & 
            obsval, obs_val(1), qc_val(1)

real(r8), allocatable ::  qc(:),obs_time(:), obs_altitude(:)
real(r8), allocatable ::  obs_longitude(:), obs_latitude(:)
real(r8), allocatable ::  wind_u(:),wind_u_oerr(:),wind_v(:),wind_v_oerr(:), &
       tn(:), tn_oerr(:), o2(:), o2_oerr(:), n2(:), n2_oerr(:), &
       o1(:), o1_oerr(:), rho(:), rho_oerr(:), ti(:), ti_oerr(:), op(:), op_oerr(:)
type(obs_def_type)      :: obs_def
type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: time_obs, prev_time, time_anal

!------------------------------------------------------------------------
!  Declare namelist parameters
!------------------------------------------------------------------------

character(len=256) :: gdc_netcdf_file     = 'gdc_input.nc'
character(len=256) :: gdc_netcdf_filelist = ''
character(len=256) :: gdc_out_file        = 'obs_seq.gdc'

namelist /convert_gdc_nml/gdc_netcdf_file,          &
                                  gdc_netcdf_filelist, gdc_out_file, datestr, obs_window

! initialize some values
obs_num   = 1
first_obs = .true.
call set_calendar_type(GREGORIAN)

!  read the necessary parameters from input.nml
call initialize_utilities()

call find_namelist_in_file("input.nml", "convert_gdc_nml", iunit)
read(iunit, nml = convert_gdc_nml, iostat = io)
call check_namelist_read(iunit, io, "convert_gdc_nml")

! Record the namelist values used for the run 
if (do_nml_file()) write(nmlfileunit, nml=convert_gdc_nml)
if (do_nml_term()) write(     *     , nml=convert_gdc_nml)

! cannot have both a single filename and a list; the namelist must
! shut one off.
if (gdc_netcdf_file /= '' .and. gdc_netcdf_filelist /= '') then
  call error_handler(E_ERR, 'convert_gdc_cdf',                     &
                     'One of gdc_netcdf_file or filelist must be NULL', &
                     source, revision, revdate)
endif
if (gdc_netcdf_filelist /= '') from_list = .true.

! need to know a reasonable max number of obs that could be added here.
if (from_list) then
   call find_textfile_dims(gdc_netcdf_filelist, nfiles, dummy)
   num_new_obs = 100000 * nfiles
else
   num_new_obs = 100000
endif

!  either read existing obs_seq or create a new one
call static_init_obs_sequence()
call init_obs(obs, num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

read(datestr(1:4),   fmt='(i4)') iyear
read(datestr(6:7),   fmt='(i2)') imonth
read(datestr(9:10),  fmt='(i2)') iday
read(datestr(12:13), fmt='(i2)') ihour
read(datestr(15:16), fmt='(i2)') imin
read(datestr(18:19), fmt='(i2)') isec
time_anal = set_date(iyear, imonth, iday, ihour, imin, isec)
call get_time(time_anal, asec, aday)



inquire(file=gdc_out_file, exist=file_exist)
if ( file_exist ) then

   write(msgstring, *) "found existing obs_seq file, appending to ", trim(gdc_out_file)
   write(msgstring2, *) "adding up to a maximum of ", num_new_obs, " new observations"
   call error_handler(E_MSG, 'convert_gdc_cdf', msgstring, &
                      source, revision, revdate, text2=msgstring2)
   call read_obs_seq(gdc_out_file, 0, 0, num_new_obs, obs_seq)

else

   write(msgstring,  *) "no previously existing obs_seq file, creating ", trim(gdc_out_file)
   write(msgstring2, *) "with up to a maximum of ", num_new_obs, " observations"
   call error_handler(E_MSG, 'convert_gdc_cdf', msgstring, &
                      source, revision, revdate, text2=msgstring2)

   call init_obs_sequence(obs_seq, num_copies, num_qc, num_new_obs)
   call set_copy_meta_data(obs_seq, 1, 'observation')
   call set_qc_meta_data(obs_seq, 1, 'Data QC')

end if


! main loop that does either a single file or a list of files.
! the data is currently distributed as a single profile per file.

filenum = 1
numrejected = 0


time_obs = set_date(1970, 1, 1, 0, 0, 0)
call get_time(time_obs,  st_sec, st_day)

fileloop: do      ! until out of files

   ! get the single name, or the next name from a list
   if (from_list) then 
      next_infile = get_next_filename(gdc_netcdf_filelist, filenum)
   else
      next_infile = gdc_netcdf_file
      if (filenum > 1) next_infile = ''
   endif
   if (next_infile == '') exit fileloop
  
   call nc_check( nf90_open(next_infile, nf90_nowrite, ncid), 'file open', next_infile)
   
   call nc_check( nf90_inq_dimid(ncid,"n_timestamps", varid), 'inq dimid n_timestamps', next_infile)
   call nc_check( nf90_inquire_dimension(ncid, varid, len = timestamps),'inq dim n_timestamps',   next_infile)
   allocate(obs_time(timestamps))
   allocate(obs_longitude( timestamps))
   allocate(obs_latitude( timestamps))
   allocate(obs_altitude( timestamps))

   call nc_check( nf90_inq_varid(ncid,"geod_alt", varid)    ,'inq varid alt', next_infile)
   call nc_check( nf90_get_var(ncid, varid, obs_altitude)   ,'get var   alt', next_infile)
   call nc_check( nf90_inq_varid(ncid,"geod_lon", varid)    ,'inq varid lon', next_infile)
   call nc_check( nf90_get_var(ncid, varid, obs_longitude)  ,'get var   lon', next_infile)
   call nc_check( nf90_inq_varid(ncid,"geod_lat", varid)    ,'inq varid lat', next_infile)
   call nc_check( nf90_get_var(ncid, varid, obs_latitude)   ,'get var   lat', next_infile)
   call nc_check( nf90_inq_varid(ncid,"time", varid)        ,'inq varid time', next_infile)
   call nc_check( nf90_get_var(ncid, varid, obs_time),'get var   time', next_infile)

   allocate(wind_u( timestamps))
   allocate(wind_u_oerr( timestamps))
   allocate(wind_v( timestamps))
   allocate(wind_v_oerr( timestamps))
   allocate(tn( timestamps))
   allocate(tn_oerr( timestamps))
   allocate(o2( timestamps))
   allocate(o2_oerr( timestamps))
   allocate(n2( timestamps))
   allocate(n2_oerr( timestamps))
   allocate(o1( timestamps))
   allocate(o1_oerr( timestamps))
   allocate(rho( timestamps))
   allocate(rho_oerr( timestamps))
   
   allocate(ti( timestamps))
   allocate(ti_oerr( timestamps))
   allocate(op( timestamps))
   allocate(op_oerr( timestamps))
 
   allocate(qc( timestamps))

   call nc_check( nf90_inq_varid(ncid,"TN", varid)          ,'inq varid variable', next_infile)
   call nc_check( nf90_get_var(ncid, varid, tn)             ,'get var variable', next_infile)
   tn_oerr = tn * 0.02;
   call nc_check( nf90_inq_varid(ncid,"UN", varid)          ,'inq varid variable', next_infile)
   call nc_check( nf90_get_var(ncid, varid, wind_u)         ,'get var variable', next_infile)
   wind_u_oerr = 4.5;
   call nc_check( nf90_inq_varid(ncid,"VN", varid)          ,'inq varid variable', next_infile)
   call nc_check( nf90_get_var(ncid, varid, wind_v)         ,'get var variable', next_infile)
   wind_v_oerr = 4.5;
   call nc_check( nf90_inq_varid(ncid,"O1", varid)          ,'inq varid variable', next_infile)
   call nc_check( nf90_get_var(ncid, varid, o1)              ,'get var variable', next_infile)
   o1_oerr = o1*0.01;
   call nc_check( nf90_inq_varid(ncid,"O2", varid)          ,'inq varid variable', next_infile)
   call nc_check( nf90_get_var(ncid, varid, o2)             ,'get var variable', next_infile)
   o2_oerr = o2*0.01;
   call nc_check( nf90_inq_varid(ncid,"N2", varid)          ,'inq varid variable', next_infile)
   call nc_check( nf90_get_var(ncid, varid, n2)             ,'get var variable', next_infile)
   n2_oerr = n2*0.01;
   call nc_check( nf90_inq_varid(ncid,"Rho", varid)          ,'inq varid variable', next_infile)
   call nc_check( nf90_get_var(ncid, varid, rho)             ,'get var variable', next_infile)
   rho_oerr = rho*0.01;
   call nc_check( nf90_inq_varid(ncid,"OP", varid)          ,'inq varid variable', next_infile)
   call nc_check( nf90_get_var(ncid, varid, op)             ,'get var variable', next_infile)
   op_oerr = op*0.01;
   call nc_check( nf90_inq_varid(ncid,"TI", varid)          ,'inq varid variable', next_infile)
   call nc_check( nf90_get_var(ncid, varid, ti)             ,'get var variable', next_infile)
   ti_oerr = ti*0.02;
 
   qc = 0

 
! read variable from netcdf  need to be fixed later
!   call nc_check( nf90_inq_varid(ncid,"variable name", varid)          ,'inq varid variable', next_infile)
!   call nc_check( nf90_get_var(ncid, varid, variable)                  ,'get var variable', next_infile)
!   call nc_check( nf90_inq_varid(ncid,"observation error name", varid) ,'inq varid observation error', next_infile)
!   call nc_check( nf90_get_var(ncid, varid, variable_oerr)             ,'get var   U Error', next_infile)
!   call nc_check( nf90_inq_varid(ncid,"data quality", varid)           ,'inq varid QC', next_infile)
!   call nc_check( nf90_get_var(ncid, varid, qc)                        ,'get var   QC', next_infile)
  
   call nc_check( nf90_close(ncid) , 'close file', next_infile)
   
    print*,timestamps


   obsloop2: do k = 1, timestamps
!-----------------------------------------------------------------------------------------------------------------------------------

     qc_val(1)  = qc(k)
     if (qc_val(1) .eq. 0) then
        oday=int((obs_time(k)+st_sec)/86400)+st_day;
        osec=int(obs_time(k)+st_sec-(oday-st_day)*86400);

        time_obs = set_time(osec,oday)
        if  ( abs(real(osec+oday*86400)-(asec+aday*86400)) < obs_window/2 ) then
        ! lon(zloc) and lon(zloc+1) range from -180 to +180
        lono=obs_longitude(k)
        lato=obs_latitude(k)
        hghto=obs_altitude(K)*1000 ! km ==> m
        if (lono < 0) lono=lono+360
        call set_obs_def_location(obs_def,set_location(lono,lato,hghto,VERTISHEIGHT))
        call set_obs_def_time(obs_def, time_obs)
     
        oerr= wind_u_oerr(k)
        obs_val(1) = wind_u(k)
        call set_obs_def_type_of_obs(obs_def, GDC_VELOCITY_U)
        call set_obs_def_error_variance(obs_def, oerr * oerr)
        call set_obs_def_key(obs_def, obs_num)
        call set_obs_def(obs, obs_def)
        call set_obs_values(obs, obs_val)
        call set_qc(obs, qc_val)

        call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
        obs_num = obs_num+1

        oerr= wind_v_oerr(k)
        obs_val(1) = wind_v(k)
        call set_obs_def_type_of_obs(obs_def, GDC_VELOCITY_V)
        call set_obs_def_error_variance(obs_def, oerr * oerr)
        call set_obs_def_key(obs_def, obs_num)
        call set_obs_def(obs, obs_def)
        call set_obs_values(obs, obs_val)
        call set_qc(obs, qc_val)

        call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
        obs_num = obs_num+1

        oerr= tn_oerr(k)
        obs_val(1) = tn(k)
        call set_obs_def_type_of_obs(obs_def, GDC_TEMPERATURE)
        call set_obs_def_error_variance(obs_def, oerr * oerr)
        call set_obs_def_key(obs_def, obs_num)
        call set_obs_def(obs, obs_def)
        call set_obs_values(obs, obs_val)
        call set_qc(obs, qc_val)

        call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
        obs_num = obs_num+1

        oerr= o2_oerr(k)
        obs_val(1) = o2(k)
        call set_obs_def_type_of_obs(obs_def, GDC_DENSITY_NEUTRAL_O2)
        call set_obs_def_error_variance(obs_def, oerr * oerr)
        call set_obs_def_key(obs_def, obs_num)
        call set_obs_def(obs, obs_def)
        call set_obs_values(obs, obs_val)
        call set_qc(obs, qc_val)

        call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
        obs_num = obs_num+1

        oerr= n2_oerr(k)
        obs_val(1) = n2(k)
        call set_obs_def_type_of_obs(obs_def, GDC_DENSITY_NEUTRAL_N2)
        call set_obs_def_error_variance(obs_def, oerr * oerr)
        call set_obs_def_key(obs_def, obs_num)
        call set_obs_def(obs, obs_def)
        call set_obs_values(obs, obs_val)
        call set_qc(obs, qc_val)

        call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
        obs_num = obs_num+1

        oerr= o1_oerr(k)
        obs_val(1) = o1(k)
        call set_obs_def_type_of_obs(obs_def, GDC_DENSITY_NEUTRAL_O)
        call set_obs_def_error_variance(obs_def, oerr * oerr)
        call set_obs_def_key(obs_def, obs_num)
        call set_obs_def(obs, obs_def)
        call set_obs_values(obs, obs_val)
        call set_qc(obs, qc_val)

        call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
        obs_num = obs_num+1

        oerr= ti_oerr(k)
        obs_val(1) = ti(k)
        call set_obs_def_type_of_obs(obs_def, GDC_TEMPERATURE_ION)
        call set_obs_def_error_variance(obs_def, oerr * oerr)
        call set_obs_def_key(obs_def, obs_num)
        call set_obs_def(obs, obs_def)
        call set_obs_values(obs, obs_val)
        call set_qc(obs, qc_val)

        call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
        obs_num = obs_num+1

        oerr= op_oerr(k)
        obs_val(1) = op(k)
        call set_obs_def_type_of_obs(obs_def, GDC_DENSITY_ION_OP)
        call set_obs_def_error_variance(obs_def, oerr * oerr)
        call set_obs_def_key(obs_def, obs_num)
        call set_obs_def(obs, obs_def)
        call set_obs_values(obs, obs_val)
        call set_qc(obs, qc_val)

        call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
        obs_num = obs_num+1

        oerr= rho_oerr(k)
        obs_val(1) = rho(k)
        call set_obs_def_type_of_obs(obs_def, GDC_DENSITY)
        call set_obs_def_error_variance(obs_def, oerr * oerr)
        call set_obs_def_key(obs_def, obs_num)
        call set_obs_def(obs, obs_def)
        call set_obs_values(obs, obs_val)
        call set_qc(obs, qc_val)

        call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
        obs_num = obs_num+1

        endif
     endif
   end do obsloop2

  ! clean up and loop if there is another input file
   deallocate(obs_time,qc)
   deallocate(obs_altitude, obs_longitude, obs_latitude)
   deallocate(wind_u, wind_u_oerr, wind_v, wind_v_oerr)
   deallocate(o2, o2_oerr, n2, n2_oerr, o1, o1_oerr, rho, rho_oerr)
   deallocate(tn, tn_oerr, ti, ti_oerr, op, op_oerr)

  filenum = filenum + 1

end do fileloop

! done with main loop.  if we added any new obs to the sequence, write it out.
print*,gdc_out_file
if (obs_num > 1) call write_obs_seq(obs_seq, gdc_out_file)

! cleanup memory
call destroy_obs_sequence(obs_seq)
call destroy_obs(obs)   ! do not destroy prev_obs, which is same as obs

write(msgstring, *) 'processed ', filenum-1, ' total profiles'
call error_handler(E_MSG, 'convert_gdc_obs', msgstring, source, revision, revdate)

if (numrejected > 0) then
   write(msgstring,  *) numrejected, ' profiles rejected for bad incoming quality control'
   call error_handler(E_MSG, 'convert_gdc_obs', msgstring, source, revision, revdate)
endif

call finalize_utilities()

end program


