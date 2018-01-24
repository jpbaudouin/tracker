!#######################################################################
!################## DEFINITION OF THE MAIN SUBROUTINE ##################
!#######################################################################


!25-10-2017 Creation
!12-12-2017 1st Commit -> up to itsloop
!24-01-2018 2nd commit -> working version



! call read_tcv_card not writen
! Add within the R code the possibility to compute vorticity from winds
! Define ts_tim_0 when reading old output for previously known object


! CONTAINS :
!   The modules main with its subroutine trackmain



!#######################################################################
!#######################################################################


module main 

  use namelist_trk 
  use data_param
  use read_datafile
  use time_step
  use write_read_output
  use grid_projection
  use core_tracker
  
 
  implicit none
   
  !variables
   
  !private
  !public trackmain
   
  contains 
  
  !---------------------------------------------------------------------  
  subroutine trackmain
  ! ABSTRACT: This is the Main subroutine which will be called by R

    implicit none
    
    integer ilast_object ! index of the last storm being detected
    integer its ! the current time being treated
    integer nobj, iobj, iobj2 ! number of objects being tracked and index
    real, allocatable :: minmax_detec2(:)
    integer, allocatable :: ind(:), sorted_list(:)
    integer i, iint, ptc_i, ptc_j
    real dt, geslon, geslat, clon, clat, cval, fixlon, fixlat
    real contour_min, contour_max, rcontour_max
    logical found
    integer err_a1, err_a2, err_a3, err_a4, err_a5, err_a6, err_a7
    integer err_a8, err_a9, err_a10, err_a11, err_a12, err_a13, err_a14
    integer err_a15, err_a16, err_a17, err_a18, err_a19, err_a20
    integer err_a21, ret
     
    logical, allocatable :: mask_object_max(:,:), mask_object_min(:,:)
    logical, allocatable :: mask_object(:,:), mask_detec(:,:)
    real, allocatable :: vals_int(:), vals_fc(:)
    
                
    if (verb .ge. 2) then
      print *, 'Loading parameters and datafile'
      print *, '###############################'
      print *, ' '
    endif
    
    ! The namelist
    call read_nlists
    if (verb .ge. 2) then
      call print_namelist
    endif
    
    ! The time steps
    call readts

    ! Loading datafile
    call open_file
    call getgridinfo
    
   
  !----------------------------------------------------------------
  !-- We allocate the different arrays (Rflag%, data and object) --
  !----------------------------------------------------------------
    allocate (Rflag%fixcenter(numfield%fixcenter),stat=err_a1)
    allocate (Rflag%advection(numfield%advection),stat=err_a2)
    allocate (Rflag%intensity(numfield%intensity),stat=err_a3)
    
    allocate (data_detection(imax,jmax),stat=err_a4)
    allocate (data_fixcenter(imax,jmax,numfield%fixcenter),stat=err_a5)
    allocate (data_intensity(imax,jmax,numfield%intensity),stat=err_a6)   
    
    allocate (data_u(imax, jmax, numfield%advection), stat=err_a7)
    allocate (data_v(imax, jmax, numfield%advection), stat=err_a8)
    
    allocate (masked_out(imax, jmax), stat=err_a9)
    allocate (mask_object_max(imax, jmax), stat=err_a10)
    allocate (mask_object_min(imax, jmax), stat=err_a11)
    allocate (mask_object(imax, jmax), stat=err_a12)
    allocate (mask_detec(imax, jmax), stat=err_a13)


    allocate (prevlon(max_object), stat=err_a14)
    allocate (prevlat(max_object), stat=err_a15)
    allocate (minmax_detec(max_object), stat=err_a16)
    allocate (speed(max_object), stat=err_a17)
    allocate (direction(max_object), stat=err_a18)
    allocate (names_obj(max_object), stat=err_a19)
    
    allocate (vals_fc(numfield%intensity), stat=err_a20)
    allocate (vals_int(numfield%intensity), stat=err_a21)
    
    prevlon       = -9999.99
    prevlat       = -9999.99
    minmax_detec  = -9999.99
    speed         = -9999.99
    direction     = -9999.99
    names_obj     = ""
    beingTracked  = .false.

 
    if (err_a1 /= 0 .or. err_a2 /= 0 .or. err_a3 /= 0 .or. &
        err_a4 /= 0 .or. err_a5 /= 0 .or. err_a6 /= 0 .or. &
        err_a5 /= 0 .or. err_a6 /= 0 .or. err_a7 /= 0 .or. &
        err_a8 /= 0 .or. err_a9 /= 0 .or. err_a10 /= 0 .or. &
        err_a11 /= 0 .or. err_a12 /= 0 .or. err_a13 /= 0 .or. &
        err_a14 /= 0 .or. err_a15 /= 0 .or. err_a16 /= 0 .or. &
        err_a17 /= 0 .or. err_a18 /= 0 .or. err_a19 /= 0 .or. &
        err_a20 /= 0 .or. err_a21 /= 0) then 
      if (verb .ge. 1) then
        print *, ""
        print *, "!!! ERROR 11: in getgridinfo_netcdf"
        print *, "!!! allocating Rflag or data arrays failed"
      endif 
      STOP 11
    endif


  !------------------------------------    
  !-- The previously tracked objects --
  !------------------------------------
    call read_input(ilast_object)
    
    if (ilast_object > 0) then
      beingTracked(1:ilast_object) = .true.
    endif


  !-------------------------------------------------------------
  !-- start of the loop over all the timestep of the datafile --
  !-------------------------------------------------------------
    
    itsloop: do its = 1,nts
    
      if ( verb .ge. 2 ) then
        print *,' '
        print *,' '
        print *,'*-------------------------------------------*'
        print *,'*-------------------------------------------*'
        print "(A, A)", "    New forecast hour: ", print_ts_tim(its)
!        call date_and_time (big_ben(1),big_ben(2),big_ben(3)
!     &                      ,date_time)
!        write (6,31) date_time(5),date_time(6),date_time(7)
! 31     format (1x,'TIMING: before stormloop  ',i2.2,':',i2.2,':',i2.2)
        print *,'*-------------------------------------------*'
        print *,'*-------------------------------------------*'
      endif
    
      ! initialisation of variable arrays and Rflag
      data_detection = -9999.0
      data_fixcenter = -9999.0
      data_intensity = -9999.0
      data_u         = -9999.0
      data_v         = -9999.0
      
      Rflag%fixcenter = .false.
      Rflag%advection = .false.
      Rflag%intensity = .false.
      
      masked_out = .false.
      mask_detec = .false.
    
    
      call getdata(ts_tim(its))
      
    !------------------------------------
    !-- sorting object by minmax_detec --
    !------------------------------------
      
      ! This part is made in order to reduce minmax_detec size to the
      !   actual nunber of object being track
      
      nobj = count(beingTracked)
      if (nobj == 0) then
        go to 9 ! no object detected, go to first_ges_center
      endif
      
      if (allocated(minmax_detec2)) deallocate (minmax_detec2)
      if (allocated(sorted_list)) deallocate (sorted_list)
      if (allocated(ind)) deallocate (ind)
      allocate(minmax_detec2(nobj))
      allocate(sorted_list(nobj))
      allocate(ind(nobj))
      
      i = 0
      do iobj=1,ilast_object
        if(beingTracked(iobj)) then
          i = i +1
          minmax_detec2(i) = minmax_detec(iobj)
          ind(i) = iobj
        endif
      enddo
      
      call qsort(minmax_detec2, sorted_list, nobj)
      
      if (fname%detec_minmax == "max") then
        sorted_list = sorted_list(nobj:1:-1)
      endif

    !------------------------------------------------------------
    !-- start of the loop over all the object already detected --
    !------------------------------------------------------------
      object_loop: do iobj2 = 1,nobj !loop on all the object detected so far
        iobj = ind(sorted_list(iobj2))
      
        if (beingTracked(iobj)) then

          ! Compute the projection of the next point from prevlon, 
          !   prevlat, speed and direction
          if (its == 1) then
            dt = ts_tim(1) - ts_tim_0
          else
            dt = ts_tim(its) - ts_tim(its-1)
          endif

          ! to get dt in seconds
          if (fileinfo%time_unit == 'minutes') then
            dt = dt*60
          else if (fileinfo%time_unit == 'hours') then
            dt = dt*3600
          else if (fileinfo%time_unit == 'days') then
            dt = dt*3600*24
          else if (fileinfo%time_unit == 'months') then
            dt = dt*3600*24*30.4375
          else if (fileinfo%time_unit == 'years') then
            dt = dt*3600*24*365.25
          endif
          

          !Note: speed is in m/s, and speed*dt/1000 is in km
          call calc_next_point(prevlon(iobj), prevlat(iobj), &
                  speed(iobj)*dt/1000, direction(iobj), geslon, geslat)
          if ( verb .ge. 3 ) then
            print *,' '
            print *,' '
            print *, "*-------------------------------------*"
            print *, "A new predicted point has been computed"
            print *, "Object:    ", names_obj(iobj)
            print *, "Time step: ", print_ts_tim(its)
            print *, "lon:       ", geslon
            print *, "lat:       ", geslat
            print *, "*-------------------------------------*"
          endif
          
          
          ! Look for a minmax in the detection field
          call find_minmax(geslon, geslat, data_detection, &
                fname%detec_minmax, ptc_i, ptc_j, clon, clat, cval, &
                found)

          if (.not. found) then
            if ( verb .ge. 3 ) then
              print *,' '
              print *, "When looking for the next position of ", &
                          names_obj(iobj)
              print *, "No minmax was found nearby the predicted point"
              print *, "The tracking of the object stops"
            endif
            beingTracked(iobj) = .false.
            cycle
          endif
                

          ! Look for a closed contour                
          if (fname%detec_minmax == 'min') then
            contour_min = cval + trkrinfo%contint
          else
            contour_min = cval - trkrinfo%contint
          endif
    
          call check_contour(ptc_i, ptc_j, data_detection, &
                fname%detec_minmax, contour_min, .false., mask_object_max, &
                mask_object_min, contour_max, rcontour_max, ret)
          
          if (ret /= 0) then
            if ( verb .ge. 3 ) then
              print *,' '
              print *, "When looking for the next position of ", &
                          names_obj(iobj)
              print *, "The new minmax found doesn't have a "// &
                          "minimum closed contour"
              print *, "The tracking of the object stops"
            endif
            beingTracked(iobj) = .false.
            cycle
          endif

          if (verb >= 3) then
            print *, "A closed contour was found"
            print *, "contour_max: ", contour_max
            print *, "rcontour_max: ", rcontour_max
          endif
          
       
              
          ! better approximation of the center                  
          call fixcenter(ptc_i, ptc_j, mask_object_max, &
                mask_object_min, mask_object, fixlon, fixlat, &
                vals_fc, contour_max, rcontour_max)


          ! compute the minmax for the fields intensity
          if (verb >= 3) then
            print *, ""
            print *, "Intensity minmax:"
            print *, "-----------------"
          endif
          do iint = 1, numfield%intensity
            if (Rflag%intensity(iint)) then
              call find_minmax(fixlon, fixlat, &
                      data_intensity(:,:,iint),fname%int_minmax(iint), &
                      ptc_i, ptc_j, clon, clat, vals_int(iint), found)
            else
              vals_int(iint) = -9999.
            endif
          enddo

          ! compute speed and direction from wind speed
          call get_speed_dir(iobj, its, fixlon, fixlat, &
                               mask_object, .false.)

          prevlon(iobj) = fixlon
          prevlat(iobj) = fixlat
          minmax_detec(iobj) = cval

          ! write ouputs
          call default_output(iobj, its, contour_max, rcontour_max, &
                              vals_fc, vals_int)


          ! fill masked_out with new trues
          masked_out = masked_out .or. mask_object
          mask_detec = mask_detec .or. mask_object_max


        endif
      enddo object_loop
      
      
      9 continue ! if no object being tracked
      
      
    !--------------------------
    !-- Look for new objects --
    !--------------------------      
      if (trkrinfo%genesis_flag) then
        call first_ges_center(its, ilast_object, mask_detec)
      endif

    enddo itsloop

  end subroutine trackmain
  
end module




