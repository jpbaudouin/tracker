!#######################################################################
!##### DEFINITION OF THE SUBROUTINES TO READ AND WRITE THE OUTPUT ######
!#######################################################################


!12-01-2018 Creation and 1st commit
!15-01-2018 2nd commit -> adding speed and direction in output
!                      -> removing file_default_output as a module variable
!24-01-2018 2nd commit -> working version


!Note: speed and direction are useful when a new instance of the program
!   is launched with previously tracked object

! CONTAINS :
!   The module "write_read_output". 
!   It includes variables saving some values related to object
!   It includes the subroutines:
!       - current_datetime(its, date, time)
!       - newobject_name(its, iobj, lon, lat) 
!       - default_output(iobj, its, contour_max, rcontour_max,
!                         minmaxs_fixcenter, minmaxs_intensity)
!       - read_input(iobj)
!       And the function "print_ts_tim(its)"




!#######################################################################
!######## Module for the reading and the writing of the output #########
!#######################################################################
!-----------------------------------------------------------------------

module write_read_output


  use namelist_trk
  use time_step


  implicit none
  
  ! vectors saving some object values. Put here necause needed in output
  !   printing, but mainly used in core module
  character*29, save, allocatable :: names_obj(:)
  real, save, allocatable :: prevlon(:), prevlat(:) !last position of any given objects 
  real, save, allocatable :: speed(:), direction(:) !projection for that last position
  real, save, allocatable :: minmax_detec(:) ! minmax in detection field at the last position
  
  integer ts_tim_0 ! from the reading of the input

  contains


  !---------------------------------------------------------------------------------
  !------- subroutine that compute date and time from startdate and timestep -------
  !---------------------------------------------------------------------------------
  
  subroutine current_datetime(its, date, time) 
  ! From trkrinfo%startdate, this subroutine adds the amount of time 
  !   given by ts_time(its) to compute the date and time associated with
  !   the timestep its.
  
  !  Note: ts_tim > 0
  
    ! INPUT
    ! its     the timestep considered
  
    ! OUTPUT
    ! date    the date of the timestep: yyyymmdd
    ! time    the time of the timestep: hhmmss
  
    implicit none
  
    integer, intent(in) :: its
    character*8, intent(out) :: date
    character*6, intent(out) :: time
  
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
  
    if (s >= 60) then
      q = s/60
      s = s - q*60
      mn = mn + q
    endif
      
    if (mn >= 60) then
      q = mn/60
      mn = mn - q*60
      h = h + q
    endif  
  
    if (h >= 24) then
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
          print *, "!!! ERROR 241 in object_name:"
          print *, "!!! the month in trkrinfo%startdate is < 0"
        endif  
        STOP 241
      endif
    enddo          
          
          
    if (m > 12) then
      q = m/12
      m = m - q*12
      y = y + q
    endif
    
    
    
    ! write the outputs
    write(date, "(I0.4, 2I0.2)")  y, m, d
    write(time, "(3I0.2)")  h, mn, s
    
    
  end subroutine
    
    
  
  
  !---------------------------------------------------------------------------------
  !-------- function outpuring a formated date and time string of character --------
  !---------------------------------------------------------------------------------
  function print_ts_tim(its)
  ! Using the subroutine current_datetime, it returns a string of
  !   character for printing date and time.
  ! The format is yyyy-mm-dd hh:mm:ss, the same as in the namelist
  
  
    implicit none
    
    integer, intent(in) :: its
    character*19 print_ts_tim
    
    character*8 date
    character*6 time
    
    call current_datetime(its, date, time)
    
    write(print_ts_tim, "(11A)") date(1:4), "-", date(5:6), "-", date(7:8), &
                        " ", time(1:2), ":", time(3:4), ":", time(5:6)
                        
  end function
    
  
  
    
  
  !---------------------------------------------------------------------------------
  !-------- subroutine that fill names_obj with the name of a new object -----------
  !---------------------------------------------------------------------------------
  
  subroutine newobject_name(its, iobj, lon, lat) 
  ! Each object have an identifier (an integer) to which is associate a
  !   name that will be print out in the output as well as in the user
  !   messages. The name has the format yyyymmddhhmmss_lonc_latc, where
  !   yyyymmddhhmmss is the start date and time of the object track, and
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
    
    character*8  date
    character*6  time
    character    lonc, latc
    integer         lat2, lon2
    character*29 object_name
    
  
  
  !------------------------------------------------------------
  !-- calculate the starting date and time of the new object --
  !------------------------------------------------------------
    call current_datetime(its, date, time) 

  !--------------------------
  !-- check lon/lat format --
  !--------------------------
  
    if (lon >= 180 .and. lon < 360) then
      lon2 = (360-lon)*100
      lonc = "W"
    else if (lon < 0 .and. lon >= -180) then
      lon2 = -lon*100
      lonc = "W"
    else if (lon >= 0 .and. lon < 180) then
      lon2 = lon*100
      lonc = "E"
    else
      if ( verb .ge. 1 ) then
        print *, " "
        print *, "!!! ERROR 251 in newobject_name:"
        print *, "!!! The longitude exceeded the limit [-180,260]"
        print *, "!!! time step: ", its
        print *, "!!! object:    ", iobj
        print *, "!!! lon:       ", lon
        print *, "!!! lat:       ", lat
      endif  
      STOP 251
    endif

    
    if (lat < 0 .and. lat >= -90) then
      lat2 = -lat*100
      latc = "S"
    else if (lat >= 0 .and. lat <= 90) then
      lat2 = lat*100
      latc = "N"
    else
      if ( verb .ge. 1 ) then
        print *, " "
        print *, "!!! ERROR 252 in newobject_name:"
        print *, "!!! The latitude exceeded the limit [-90,90]"
        print *, "!!! time step: ", its
        print *, "!!! object:    ", iobj
        print *, "!!! lon:       ", lon
        print *, "!!! lat:       ", lat
      endif  
      STOP 252
    endif
  
    
  !---------------------------
  !-- create the new string --
  !---------------------------
    write(object_name, "(A8, A, A6, A, I0.5, A, A, I0.5, A)") date, &
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
  subroutine default_output(iobj, its, contour_max, rcontour_max, &
                    minmaxs_fixcenter, minmaxs_intensity)
                 
  ! Note: iobj give the object's name through names_obj               
                 
                    
    implicit none
    
    integer, intent(in) :: iobj, its
    real,    intent(in) :: minmaxs_fixcenter(:), minmaxs_intensity(:)
    real,    intent(in) :: contour_max, rcontour_max
  
    integer i
    character*8 date
    character*6 time
    character*255 header, obj_info, file_default_output, file_last_ts
    logical OK


    if ( verb .ge. 3 ) then
      print *, ""
      print *, "printing the outputs ..."
    endif
    
    file_default_output = "output_"// trim(fileinfo%model_name)//"_"// &
                            trim(trkrinfo%run_name)//".txt"

    
    inquire(unit=64, opened=OK)
    if (.not. OK) then
      open(file=trim(file_default_output), unit = 64)
      
      print *, "print header"
      
      header = "  ID,                          name, yyyymmdd, "// &
               "hhmmss,    lon,    lat,       speed,  dir, "

      write(header, "(A, A12, A)") trim(header), &
                                  adjustr(trim(fname%detection)), &
                                  ", contour_max,    rad_max, "
                                  




      do i = 1,numfield%fixcenter
        write(header, "(A, A11, A)") trim(header), &
                                    adjustr(trim(fname%fixcenter(i))), ", "
      enddo

      do i = 1,numfield%intensity
        write(header, "(A, A11, A)") trim(header), &
                                    adjustr(trim(fname%intensity(i))), ", "
      enddo
      
      write(64, "(A)") trim(header)


    endif


    call current_datetime(its, date, time)   

    
    write(obj_info, "(I4, A, A29, A, A8, A, A6, A, F6.2, A, F6.2, " // &
                      "A, EN10.2, A, F4.0, A, ES11.4, A, ES11.4, " // &
                      "A, EN10.2, A)") &
                  iobj, ", ", names_obj(iobj), ", ", date, ", ", &
                  time, ", ", prevlon(iobj), ", ", prevlat(iobj), &
                  ", ", speed(iobj), ", ", direction(iobj), ", ", &
                  minmax_detec(iobj), ", ", contour_max, ", ", &
                  rcontour_max, ", "
      
    do i = 1,numfield%fixcenter
      write(obj_info, "(A, ES11.4, A)") trim(obj_info), &
                                        minmaxs_fixcenter(i), ", "
    enddo
      
    do i = 1,numfield%intensity
      write(obj_info, "(A, ES11.4, A)") trim(obj_info), &
                                        minmaxs_intensity(i), ", "
    enddo
    
    write(64, "(A)") trim(obj_info)


    flush(64)
    
    
    ! At the last time step, the objects info are also printed in a
    !   separate files to be potentially read by read_output in an other
    !   instance of the programme (mainly of use when datafiles are split)
    if (its == nts) then

      file_last_ts = "output_"// trim(fileinfo%model_name)//"_"// &
                        trim(trkrinfo%run_name)//"_last_ts.txt"
                            
      inquire(unit=65, opened=OK )
      if (.not. OK) then
        open(file=trim(file_last_ts), unit = 65)
      endif

      write(obj_info, "(I4, A, A29, A, A8, A, A6, A, F6.2, A, F6.2, " // &
                      "A, EN10.2, A, F4.0, A, ES11.4)") &
                          iobj, ", ", names_obj(iobj), ", ", date, &
                          ", ", time, ", ", prevlon(iobj), ", ", &
                          prevlat(iobj),  ", ", speed(iobj), ", ", &
                          direction(iobj), ", ", minmax_detec(iobj)

      write(65, "(A)") trim(obj_info)


      flush(65)
      
    endif

    
    
    end subroutine  
  
  
  
  !---------------------------------------------------------------------------------
  !-------------------- Subroutine reading possible object inputs ---------------------
  !---------------------------------------------------------------------------------
  subroutine read_input(iobj)
  ! Read a file that was produced at the last time step of a previous 
  !   run of this program (see above default_output).
  ! This is meant to connect tracks when datafiles are split over 
  !   different time period
  
    ! OUTPUT
    ! iobj    number of object discovered in the input file


    implicit none
    
    integer, intent(out) :: iobj
    
    character*100 file_last_ts
    logical there, is_neg
    integer ierr
    integer y, m, d, h, mn, s ! date and time of the startdate
    integer y0, m0, d0, h0, mn0, s0 ! the date and time for the objects
    integer y1, m1, d1, h1, mn1, s1 ! when reading inputs and earliest date
                                    ! between datartdate and objects date
    integer y2, m2, d2 ! latest date
    integer d_m ! number of day in a given month

    integer iobj2 ! not used
    character c(5) ! not used
    


    if ( verb .ge. 3 ) then    
      print *, ""
      print *, ""
      print *, "Checking possible tracks input:"
      print *,'--------------------------------'
    endif
    
    iobj = 0
    file_last_ts = "output_"// trim(fileinfo%model_name) // "_" // &
                      trim(trkrinfo%run_name) // "_last_ts.txt"


    inquire(file=trim(file_last_ts), exist=there)
    if (.not. there) then
      if ( verb .ge. 3 ) then
        print *, ""
        print *, "There is no files with previously detected objects"
      endif
      return
    endif


    open(file=trim(file_last_ts), unit = 63)

    
    do while(.true.) ! until read reach end of file
      
      iobj = iobj + 1
      
      ! Note: the A between F7.0 and F7.2 towards the end includes all 
      !         fixcenter/intensity minmax and contour values
      
                      
      read(63, "(I4, 2X, A29, 2X, I4, I2, I2, 2X, I2, I2, I2, 2X, " // &
               "F6.2, 2X, F6.2, 2X, EN10.2, 2X, F4.0, 2X, ES11.4)", &
           iostat = ierr, end = 269)  iobj2, names_obj(iobj), &
            y1, m1, d1, h1, mn1, s1, prevlon(iobj), prevlat(iobj), &
            speed(iobj), direction(iobj), minmax_detec(iobj)
            
   !-----------------------
   !-- check values read --
   !-----------------------
      if (ierr /= 0) then
        if ( verb .ge. 1 ) then
          print *, " "
          print *, "!!! ERROR 271 in read_output:"
          print *, "!!! There was an error in the format "
          print *, "!!! when reading line ", iobj
        endif  
        STOP 271
      endif
    
      if (iobj == 1) then
        if(m1<0 .or. m1>12 .or. d1<0 .or. d1>31 .or. h1<0 .or. h1>24 &
                .or. mn1<0 .or. mn1>60 .or. s1<0 .or. s1>60) then
          if ( verb .ge. 1 ) then
            print *, " "
            print *, "!!! ERROR 272 in read_output:"
            print *, "!!! There was an error in the date format"
            print *, "!!! when reading line ", iobj
            print *, "!!! date = yyyymmdd, time = hhmmss"
            print *, "!!! year:   ", y1
            print *, "!!! month:  ", m1
            print *, "!!! day:    ", d1
            print *, "!!! hour:   ", h1
            print *, "!!! minute: ", mn1
            print *, "!!! second: ", s1
          endif  
          STOP 272
        endif
        
        y0 = y1
        m0 = m1
        d0 = d1
        h0 = h1
        mn0 = mn1
        s0 = s1
        
      else if (y0/=y1 .or. h0/=h1 .or. d0/=d1 .or. h0/=h1 .or. &
              mn0/=mn1 .or. s0/=s1) then
        if ( verb .ge. 1 ) then
          print *, " "
          print *, "!!! ERROR 273 in read_output:"
          print *, "!!! There was an error in the date"
          print *, "!!! when reading line ", iobj
          print *, "!!! All lines should have the same date and time"
        endif  
        STOP 273
      endif
      
      if (prevlon(iobj) > 360 .or. prevlon(iobj) < 0 .or. &
          prevlat(iobj) > 90 .or. prevlat(iobj) < -90) then
        if ( verb .ge. 1 ) then
          print *, " "
          print *, "!!! ERROR 274 in read_output:"
          print *, "!!! There was an error in the coordinate"
          print *, "!!! when reading line ", iobj
          print *, "!!! longitude should be between 0 and 360"
          print *, "!!! latitude should be between -90 and 90"
          print *, "!!! lon: ", prevlon(iobj)
          print *, "!!! lat: ", prevlat(iobj)
        endif  
        STOP 274
      endif
  
      if (direction(iobj) > 360 .or. direction(iobj) < 0) then
        if ( verb .ge. 1 ) then
          print *, " "
          print *, "!!! ERROR 275 in read_output:"
          print *, "!!! There was an error in the direction"
          print *, "!!! when reading line ", iobj
          print *, "!!! direction should be between 0 and 360"
          print *, "!!! direction: ", direction(iobj)
        endif  
        STOP 275
      endif


      if ( verb .ge. 3 ) then
        print *, ""
        print *, "A new object has been read in the input file: "
        print *, "name:      ", names_obj(iobj)
        print *, "longitude: ", prevlon(iobj)
        print *, "latitude:  ", prevlat(iobj)
        print *, "intensity: ", minmax_detec(iobj), " (detection field)"
        print *, "speed:     ", speed(iobj)
        print *, "direction: ", direction(iobj)
      endif

    enddo
    
    269 continue
    
    iobj = iobj - 1
    
    close(63)



                           

 !----------------------
 !-- compute ts_tim_0 --
 !----------------------

    read(trkrinfo%startdate, "(I4,A,I2,A,I2,A,I2,A,I2,A,I2)") &
      y, c(1:1), m, c(2:2), d, c(3:3), h, c(4:4), mn, c(5:5), s
  

    ts_tim_0 = y0 - y
    if (fileinfo%time_unit == 'years') then  
      return
    endif
    
    ts_tim_0 = ts_tim_0*12 + m0 - m
    if (fileinfo%time_unit == 'months') then  
      return
    endif


    if (ts_tim_0 < 0) then
      ! input file has the earliest date
      y1 = y0
      m1 = m0
      d1 = d0
      
      y2 = y
      m2 = m
      d2 = d
      
      is_neg = .true.
    
    else if (ts_tim_0 > 0) then
      ! startdate is the earliest date
      y1 = y
      m1 = m
      d1 = d
      
      y2 = y0
      m2 = m0
      d2 = d0

      is_neg = .false.
      
    else ! the two date are in the same month of the same year
    
      ts_tim_0 = d0 - d
      go to 270
      
    endif


    iobj = 0
    ! we compute the number of day between (y1, m1, d1) and (y2, m2, d2)
    ts_tim_0 = 0
    do while (y1 /= y2 .or. m1 /= m2)
      iobj = iobj + 1
      if (m1 == 1 .or. m1 == 3 .or. m1 == 5 .or. m1 == 7 .or. &
          m1 == 8 .or. m1 == 10 .or. m1 == 12) then
        d_m = 31
      else if (m1 == 4 .or. m1 == 6 .or. m1 == 9 .or. m1 == 11) then
        d_m = 30
      else if (m1 == 2) then
        if (y1 == (y1/400)*400 .or. y1 == (y1/4)*4 .and. &
            y1 /= (y1/100)*100) then
          d_m = 29
        else
          d_m = 28
        endif
      else if (m1 == 13) then
        m1 = 1
        y1 = y1 + 1
        cycle
      endif
          
      ts_tim_0 = ts_tim_0 + d_m - d1 + 1
      d1 = 1
      m1 = m1 + 1
    enddo
    
    ts_tim_0 = ts_tim_0 + d2 - d1

    if (is_neg) then
      ts_tim_0 = -ts_tim_0
    endif
    
    270 continue ! same month and year    
    if (fileinfo%time_unit == 'days') then  
      return
    endif
    
    ts_tim_0 = ts_tim_0 * 24 + h0 - h
    if (fileinfo%time_unit == 'hours') then  
      return
    endif

    ts_tim_0 = ts_tim_0 * 60 + mn0 - mn
    if (fileinfo%time_unit == 'minutes') then  
      return
    endif 

    ts_tim_0 = ts_tim_0 * 60 + s0 - s
    ! only fileinfo%time_unit = 'seconds' remains
    
    
  end subroutine
    

end module
  
  
  




