! This aims to test all the subroutine from the tracker

!Test stops

program test_Tracker

  use namelist_trk 
  use data_param
  use time_step
  use write_read_output
  use grid_projection
  use core_tracker
  use main
    
  implicit none
  
  character*7 save_time_unit
  real dist, dir, lon2, lat2
  integer iobj
  
  character*8 date
  character*6 time
  
  
  allocate(ts_tim(10))
  allocate(prevlon(10))
  allocate(prevlat(10))
  allocate(speed(10))
  allocate(direction(10))
  allocate(names_obj(10))
  allocate(intensity(10))
  
  
  
  ! namelist_trk from namelist.f90
  print *, "Test namelist reading"
  call read_nlists
  
  
  ! write_read_output from write_read_output.f90
  save_time_unit = fileinfo%time_unit

  
  ts_tim(1) = 0
  ts_tim(2) = 10
  ts_tim(3) = 100
  ts_tim(4) = 1096
  ts_tim(5) = 10000
  ts_tim(6) = 10
  ts_tim(7) = 100
  ts_tim(8) = 100

  print *, "test date/time conversion in current_datetime()"
  print *, " and printing with print_ts_tim()"
  
  print *, "Same as stardate: ", print_ts_tim(1)
  fileinfo%time_unit = "years"
  print *, "+10 years:        ", print_ts_tim(2)
  fileinfo%time_unit = "months"
  print *, "+100 months:      ", print_ts_tim(3)
  fileinfo%time_unit = "days"
  print *, "+1096 days (1997):", print_ts_tim(4)
  fileinfo%time_unit = "hours"
  print *, "+10000 hours:     ", print_ts_tim(5)
  print *, "+10 hours:        ", print_ts_tim(6)
  fileinfo%time_unit = "minutes"
  print *, "+100 minutes:     ", print_ts_tim(7)
  fileinfo%time_unit = "seconds"
  print *, "+100 seconds:     ", print_ts_tim(8)

  fileinfo%time_unit = save_time_unit
  
  call newobject_name(10, 1, 0., 0.) ! its, iobj, lon, lat
  print *, "object name at ts 10, lon 0, lat 0: ", names_obj(1)
  call newobject_name(10, 2, -10.546871687, -59.) ! its, iobj, lon, lat
  print *, "object name at ts 10, lon -10.54..., lat -59: ", names_obj(2)
  
  nts = 10
  
  prevlat(1) = 5.5
  prevlon(1) = 6.5
  speed(1)   = 1.1
  direction(1) = 180.
  intensity(1) = -15.
  
  prevlat(2) = -32.5
  prevlon(2) = 112.84563
  speed(2)   = 96.1
  direction(2) = 0.98
  intensity(2) = -0.097
  call default_output(1, 10, 100., 1000., (/-1.,-2./), (/-10./))
  call default_output(2, 10, 100., 1000., (/-78.,-2.09/), (/-110./))
  
  !call system("tail -2 output_"// trim(fileinfo%model_name) //"_"// &
  !            trim(trkrinfo%run_name)   //".txt > output_"// &
  !            trim(fileinfo%model_name) //"_"// &
  !            trim(trkrinfo%run_name)   //"_last_ts.txt")
  
  close(65)
  
  call read_input(iobj)
  print *, "number of object read in input (2)", iobj
  print *, "initial time: ", ts_tim(10), ", compare to: ", ts_tim_0
  
  ! grid_projection from grid_projection.f90
  print *, "test calc_dist_dir"
  call calc_dist_dir(0.,0.,0.,0.,.true.,dist, dir)
  print *, "same point: ", dist, dir
  call calc_dist_dir(-5.,0.,5.,0.,.true.,dist, dir)
  print *, "distance = 1111.95 (1/10 of earth circum along eq): ", dist
  print *, "direction = 90 (East, accross Greenwich)", dir
  call calc_dist_dir(185.,0.,175.,0.,.true.,dist, dir)
  print *, "Same but accross 180° and westwards (270): ", dist, dir
  call calc_dist_dir(0.,0.,0.,90.,.true.,dist, dir)
  print *, "distance = 10007.55 (from eq to pole): ", dist
  print *, "direction = 0 (North)", dir
  call calc_dist_dir(0.,90.,0.,-90.,.true.,dist, dir)
  print *, "distance = 20015.1 (from pole to pole): ", dist
  print *, "direction = 180 (South)", dir
  call calc_dist_dir( 0.1236,52.1588,4.6994,46.8583, .true.,dist, dir)
  print *, "distance = 674.6 (Cambridge - Santenay): ", dist
  print *, "direction = 152.4483 ", dir

  print *, ""
  print *, "test calc_next_point"
  call calc_next_point(0.,0.,0.,0.,lon2,lat2)
  print *, "same point (0,0): ", lon2,lat2
  call calc_next_point(-5.,0.,1111.95, 90.,lon2,lat2)
  print *, "point (5,0) along eq, accross Greenwich: ", lon2,lat2
  call calc_next_point(185.,0.,1111.95, 270.,lon2,lat2)
  print *, "Same but accross 180° and westwards (175,0): ", lon2,lat2
  call calc_next_point(0.,0.,10007.54, 0.,lon2,lat2)
  print *, "point (0,90) from eq to pole: ", lon2,lat2
  call calc_next_point(0.,0.,20015.08 , 180.,lon2,lat2)
  print *, " point (180, 0) from eq to eq through S pole: ", lon2,lat2
  call calc_next_point(4.6994,46.8583 ,674.6,152.4483-180.,lon2,lat2)
  print *, "Cambridge (0.1236,52.1588) from Santenay: ", lon2,lat2

  
  !-----------------------
  !-- Part needing data --
  !-----------------------
  
  ! netcdf_datafile from open_datafiles.f90
  ! call open_netcdf_file
  ! call getgridinfo_netcdf
  ! call getdata_netcdf(ts)
  ! call find_var(nvar, name_var, level, ts, readflag, data_var)
  
  
  ! time_step from open_datafiles.f90
  ! call readfort15
  ! call readts_netcdf

 
end program







