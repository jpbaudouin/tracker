!#######################################################################
!##### DEFINITION OF THE SUBROUTINES TO READ AND WRITE THE OUTPUT ######
!#######################################################################


!12-01-2018 Creation and 1st commit



! CONTAINS :
!   The module "read_write_output" with the subroutines:
! 




!#######################################################################
!######## Module for the reading and the writing of the output #########
!#######################################################################
!-----------------------------------------------------------------------

module read_write_output


  use namelist_trk
  use time_step


  implicit none
  
  character(len=70), save :: file_default_output   ! "output_"// trim(fileinfo%model_name)//"_"// trim(trkrinfo%run_name)//".txt"
  character*27, save, allocatable :: names_obj(:)


  contains


  !---------------------------------------------------------------------------------
  !------- subroutine that compute date and time from startdate and timestep -------
  !---------------------------------------------------------------------------------
  
  subroutine current_datetime(its, date, time) 
  ! From trkrinfo%startdate, this subroutine adds the amount of time 
  !   given by ts_time(its) to compute the date and time associated with
  !   the timestep its.
  
  ! INPUT
  ! its     the timestep considered
  
  ! OUTPUT
  ! date    the date of the timestep: yymmdd
  ! time    the time of the timestep: hhmmss
  
  implicit none
  
  integer, intent(in) :: its
  character*6, intent(out) :: date, time
  
  integer y, m, d, h, mn, s, q
  character c(5)
  
  
    read(trkrinfo%startdate, "(I4,A,I2,A,I2,A,I2,A,I2,A,I2)") y, c(1:1), &
       m, c(2:2), d, c(3:3), h, c(4:4), mn, c(5:5), s
  
  
    if (fileinfo%time_unit == 'seconds') then
      s = s + ts_tim(its)
    else if (fileinfo%time_unit == 'minutes') then
      mn = mn + ts_tim(its)
    else if (fileinfo%time_unit == 'hours') then
      h = h + ts_tim(its)
    else if (fileinfo%time_unit == 'days') then
      d = d + ts_tim(its)
    else if (fileinfo%time_unit == 'months') then
      m = m + ts_tim(its)
    else if (fileinfo%time_unit == 'years') then
      y = y + ts_tim(its)
    endif
  
    if (s > 60) then
      q = s/60
      s = s - q*60
      mn = mn + q
    endif
      
    if (mn > 60) then
      q = mn/60
      mn = mn - q*60
      h = h + q
    endif  
  
    if (h > 24) then
      q = h/24
      h = h - q*24
      d = d + q
    endif  
  
    do while (d > 28)
      if (m == 1 .or. m == 3 .or. m == 5 .or. m == 7 .or. m == 8 .or. &
          m == 10 .or. m == 12) then
        if (d <= 31) then
          exit
        else
          d = d - 31
          m = m + 1
        endif
      else if (m == 4 .or. m == 6 .or. m == 9 .or. m == 11) then
        if (d <= 30) then
          exit
        else
          d = d - 30
          m = m + 1
        endif
      else if (m == 13) then
        m = 1
        y = y + 1
      else if (m > 13) then ! if month has been change, day is OK
        exit
      else if (m == 2) then
        if (y == (y/400)*400 .or. y == (y/4)*4 .and. y /= (y/100)*100) then
          if (d == 29) then
            exit
          else
            d = d - 29
            m = m + 1
          endif 
        else
          d = d - 28
          m = m + 1
        endif
      else
        if ( verb .ge. 1 ) then
          print *, " "
          print *, "!!! ERROR 55 in object_name:"
          print *, "!!! the month in trkrinfo%startdate is < 0"
        endif  
        STOP 55
      endif
    enddo          
          
          
    if (m > 12) then
      q = m/24
      m = m - q*24
      y = y + q
    endif
    
    
    
    ! write the outputs
    write(date, "(3I2)")  y, m, d
    write(time, "(3I2)")  h, mn, s
    
    
  end subroutine
    
    
  
  
  
  !---------------------------------------------------------------------------------
  !-------- subroutine that fill names_obj with the name of a new object -----------
  !---------------------------------------------------------------------------------
  
  subroutine newobject_name(its, iobj, lon, lat) 
  ! Each object have an identifier (an integer) to which is associate a
  !   name that will be print out in the output as well as in the user
  !   messages. The name has the format yymmddhhmmss_lonc_latc, where
  !   yymmddhhmmss is the start date and time of the object track, and
  !   lonc, latc its coordinate. lon and lat are > 0, with the character
  !   "c" indicating the direction (W/E and N/S).
  ! NOTE: this function is only called in first_ges_center when a new 
  !   object is found
  
  ! INPUT
  ! its     timestep when the object was first detected
  ! iobj    identifier of the object
  ! lon     longitude coordinate of the first detection
  ! lat     latitude coordinate of the first detection
  
    implicit none
    
    integer, intent(in) :: its, iobj
    real,    intent(in) :: lon, lat
    
    character*6  date, time
    character    lonc, latc
    integer      lat2, lon2
    character*27 object_name
  
  
  !------------------------------------------------------------
  !-- calculate the starting date and time of the new object --
  !------------------------------------------------------------
    call current_datetime(its, date, time) 

  !--------------------------
  !-- check lon/lat format --
  !--------------------------
  
    if (lon > 180) then
      lon2 = 360-lon
      lonc = "W"
    else if (lon < 0) then
      lon2 = -lon
      lonc = "W"
    else
      lon2 = lon
      lonc = "E"
    endif
    
    if (lat < 0) then
      lat2 = -lat
      latc = "S"
    else
      lat2 = lat
      latc = "N"
    endif
  
    
  !---------------------------
  !-- create the new string --
  !---------------------------

    write(object_name, "(A6, A, A6, A, F5.2, A, A, F5.2, A)") date, &
              "_", time, "_", lon2, lonc, "_", lat2, latc
          
          
    names_obj(iobj) = object_name
          
    if ( verb .ge. 3 ) then
      print *, " "
      print *, " The new object's name is: "
      print *, object_name
    endif


          
          
  end subroutine
  

  
  !---------------------------------------------------------------------------------
  !-------------------- Subroutine writing the default outputs ---------------------
  !---------------------------------------------------------------------------------
  subroutine default_output(id_obj, its, lon, lat, &
                    minmax_detection, contour_max, rcontour_max, &
                    minmaxs_fixcenter, minmaxs_intensity)
                 
  ! Note: id_obj give the object's name through names_obj               
                 
                    
    implicit none
    
    integer, intent(in) :: id_obj, its
    real,    intent(in) :: lon, lat, minmax_detection
    real,    intent(in) :: minmaxs_fixcenter(:), minmaxs_intensity(:)
    real,    intent(in) :: contour_max, rcontour_max
  
    integer unit, i
    character*6 date, time
    character*200 header, obj_info


    if ( verb .ge. 3 ) then
      print *, "printing the outputs"
    endif
    
      
    inquire(file=file_default_output, number=unit)
    if (unit == -1) then
      open(file=file_default_output, unit = 64)


      header = "ID  , names                      , yymmdd, hhmmss, "// &
               "lon    , lat     , "

      write(header, "(A, A7)") trim(header), fname%detection

      do i = 1,numfield%fixcenter
        write(header, "(A, A5, A)") trim(header), fname%fixcenter(i), ", "
      enddo

      do i = 1,numfield%intensity
        write(header, "(A, A5, A)") trim(header), fname%intensity(i), ", "
      enddo

      write(64) trim(header)


    endif  
    
  
    call current_datetime(its, date, time)   

    
    write(unit) trim(obj_info)
    
    write(obj_info, "(I4, A, A27, A, I0.8, A, I0.6, A, F7.2, A, " // &
                    "F6.2, A, F7.2, A, F5.0, A, F7.0, A)") &
                      id_obj, ", ", names_obj(id_obj), ", ", date, &
                      ", ", time, ", ", lon, ", ", lat, ", ", &
                      minmax_detection, ", ", contour_max, ", ", &
                      rcontour_max, ", "
      
    do i = 1,numfield%fixcenter
      write(obj_info, "(A, f5.0, A)") trim(obj_info), &
                                        minmaxs_fixcenter(i), ", "
    enddo
      
    do i = 1,numfield%intensity
      write(obj_info, "(A, f5.0, A)") trim(obj_info), &
                                        minmaxs_intensity(i), ", "
    enddo


    flush(64)
    
    end subroutine  
  
  
end module
  
  
  




