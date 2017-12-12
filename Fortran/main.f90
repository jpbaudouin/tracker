!#######################################################################
!################## DEFINITION OF THE MAIN SUBROUTINE ##################
!#######################################################################


!25-10-2017 Creation
!12-12-2017 1st Commit -> up to itsloop


! Add within the R code the possibility to compute vorticity from winds

! CONTAINS :
!   The modules main with its subroutine trackmain



!#######################################################################
!#######################################################################


module main 

  use core
  use namelist_trk 
  use data_param; datafile; time_step
  
  implicit none
   
  !variables
   
  !private
  !public trackmain
   
  contains 
  
  !---------------------------------------------------------------------  
  subroutine trackmain
  ! ABSTRACT: This is the Main subroutine which will be called by R

    implicit none
    
     
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
         
    ! The previously tracked storms
    call read_tcv_card ! I think it should be merged with read_gen_vitals
          ! Only the naming of the storm is different, leading to one more column
          ! Add also a more simple syntax Name/lat/lon/...
          ! Merge with sort_storms_by_pressure -> ~ sort_features_by_intensity
          ! and add an option no_sort in R when user provide tracks of previous feature
    To be done coding the subroutine writing the tcv_card

    ! Loading datafile
    call open_datafile
    call getgridinfo
    
    
    !-------------------------------------------------------------------
    !-------------------------------------------------------------------
    ! We allocate the different arrays (data and Rflag%)
    !-------------------------------------------------------------------   
    allocate (Rflag%fixcenter(numfield%fixcenter),stat=err_a1)
    allocate (Rflag%advection(numfield%advection),stat=err_a2)
    allocate (Rflag%intensity(numfield%intensity),stat=err_a3)
    
    allocate (data_detection(imax,jmax),stat=err_a4)
    allocate (data_fixcenter(imax,jmax,numfield%fixcenter),stat=err_a5)
    allocate (data_intensity(imax,jmax,numfield%intensity),stat=err_a6)   
    
    allocate (data_u(imax, jmax, numfield%advection), stat=err_a7)
    allocate (data_v(imax, jmax, numfield%advection), stat=err_a8)

 
    if (err_a1 /= 0 .or. err_a2 /= 0 .or. err_a3 /= 0 .or. &
        err_a4 /= 0 .or. err_a5 /= 0 .or. err_a6 /= 0 .or. &
        err_a7 /= 0 .or. err_a8 /= 0) then 
      if (verb .ge. 1) then
        print *, ""
        print *, "!!! ERROR 11: in getgridinfo_netcdf 
        print *, "!!! allocating Rflag or data arrays failed"
      endif 
      STOP 11
    endif
    !-------------------------------------------------------------------
    !-------------------------------------------------------------------
    
    itsloop: its = 1, nts
    
      if ( verb .ge. 2 ) then
        print *,' '
        print *,'*-------------------------------------------*'
        print "(1x,A,i4)", "*   New forecast hour: ", ts_tim(its)
!        call date_and_time (big_ben(1),big_ben(2),big_ben(3)
!     &                      ,date_time)
!        write (6,31) date_time(5),date_time(6),date_time(7)
! 31     format (1x,'TIMING: before stormloop  ',i2.2,':',i2.2,':',i2.2)
        print *,'*-------------------------------------------*'
      endif
    
      ! initialisation of variable arrays and Rflag
      data_detection = -9999.0
      data_fixcenter = -9999.0
      data_intensity = -9999.0
      data_u         = -9999.0
      data_v         = -9999.0
      
      Rflag%fixcenter = .FALSE.
      Rflag%advection = .FALSE.
      Rflag%intensity = .FALSE.
    
    
      call getdata(ts_tim(its))

    !Computing Vorticity if needed
      call subtract_cor
      call rvcal
    ! 
    
    ifhloop: ifh = 1,ifhmax !loop on the time steps

      if (trkrinfo%genflag == 'y') then !change in namelist
        call first_ges_center
      endif
      
      feature_loop: sl_counter = 1,max_feature !loop on all the feature detected so far
      
        select case (stormswitch(sl_counter)) !feature status:
          ! 1: under tracking; 2: track has ended; 3: free slots
          ! maybe only one case require action
        
        case(1) !under tracking
        
          ! Compute some KR:
           call find_maxmin !for different field
           call fixcenter
           call check_closed_contour
          
          ! In the meanwhile, test if the feature is valid 
          ! if not, the tracking should be ended (stormswitch(sl_counter) = 2)
          ! if yes, continue with further analysis
          
          ! part of it depends on the feature type
          
          call output...
          
          call get_next_ges
          
        case(2) !GFDL write some output, what for?
        
        case(3) ! Nothing
        
      
      enddo feature_loop
    
    enddo ifhloop

    call baclose
    !call w3tage
  
  end subroutine trackmain
  
end module

------------------------------------------------------------------------
subroutine tracker
  !Main loops
  
 

    stormloop: sl_counter = 1,maxstorm !storms

      if (trkrinfo%type == 'midlat' .or. trkrinfo%type == 'tcgen') 
        if (sl_counter > storms_predetected )
          if (there is mslp values)
            call first_ges_center

      select case (stormswitch(ist)) !storm validity, ist = sl_counter?

      case (1)

        call check_bounds
        if (there is zeta850 values)
          call find_maxmin('zeta850')

        if (there is zeta700 values)
          call find_maxmin('zeta700')
          
        ! same for hgt850, hgt700, mslp, vorticity_surface


        if (there is u/v850 values)
          call get_uv_guess
            if (previous is success)
              call get_uv_center

        if (there is u/v700 values)
          call get_uv_guess
            if (previous is success)
              call get_uv_center

        if (there is u/v_surface values)
          call get_uv_guess
            if (previous is success)
              call get_uv_center
              
        if (stormswitch(ist) == 1) ! might have changed after find_maxmin

          call fixcenter

          if (previous is success)
            if ((trkrinfo%type == 'midlat' .or. 'tcgen'
              .and. trkrinfo%gridtype == 'regional')
              
              if (too close of the boundaries)
                if (ifh == 1) ! why the following ??
                  call output_atcfunix
                  call output_atcf_gen
                  call output_atcf_sink



          if (previous is success) ! Why a second time ??
            if (there is mslp)
              call is_it_a_storm

            if (xinp_fixlat/lon are initialised) 
              call fix_latlon_to_ij

            if (previous is success .and isastorm(1) .and.
              (trkrinfo%type == 'midlat' .or. 'tcgen' .or.
               trkrinfo%want_oci ))

              call check_closed_contour

            if (trkrinfo%type == 'tcgen' .or. 'tracker')
              if (there is V850) ! Why only V??
                call is_it_a_storm
        
          if (trkrinfo%type == 'tracker' .or. 'tcgen')
            if (isastorm(1 and 3) .and. there is Vo850 .and.
              stormswitch(ist) == 1) ! stormswitch might have chang a second time
            
              call calcdist


          if (ifh > 1 .and. stormswitch(ist) == 1)
            if (fixlon/lat are valid)
              call calcdist
        if (there is u/v_surface .and. previous is success (fixcenter) 
          .and. stormswitch(ist) == 1)
          call get_max_wind
          call getradii



        if (phaseflag == 'y' .and. stormswitch(ist) == 1)
          call get_phase

        if (structflag == 'y' .or. ikeflag == 'y')
          call get_sfc_center

        if (structflag == 'y' .and. stormswitch(ist) == 1)
          call get_wind_structure
          if (previous is sucess)
            call output_wind_structure

        if (structflag == 'y' .and. stormswitch(ist) == 1)
          call get_fract_wind_cov
          if (previous is sucess) 
            call output_fract_wind

        if (ikeflag == 'y' .and. stormswitch(ist) == 1)
          call get_ike_stats
          if (previous is sucess)
            call output_ike
            
            
        if (previous is success (fixcenter) .and. stormswitch(ist) == 1)
          call output_atcfunix
          call get_next_ges

          call get_zeta_values

          if (trkrinfo%type == 'midlat', 'tcgen' .or. 'tracker')
            call output_atcf_gen

          call output_atcf_sink

          if (inp%model == 12 .and. ifcsthour == 0)
          ! Write vitals for GFS ens control analysis
            call output_tcvitals

          call output_hfip
        
        else

          if (ifh == 1)
            call output_atcfunix

          if (trkrinfo%type == 'midlat' .or. 'tcgen')
            call output_atcf_gen


          call output_atcf_sink
          call output_hfip

          if (trkrinfo%type == 'tracker')
            call get_next_ges

        if (ifh /= ifhmax)
          call get_next_ges

          if (trkrinfo%out_vit == 'y' .and. ifh == ifhmax - 1)
            call output_gen_vitals 
 
      case (2)

        if (ifh == 1)
          call output_atcfunix
          if (trkrinfo%type == 'midlat' .or. 'tcgen')
            call output_atcf_gen

          call output_atcf_sink
          call output_hfip

      case (3)
      
 
    enddo stormloop

    if(use_per_fcst_command=='y')
c      User wants us to run a command per forecast time

      call argreplace
      call argreplace

      call run_command(trim(pfc_final),pfcret)

    ifh = ifh + 1

    if (ifh > ifhmax) exit ifhloop
    
    if (inp%file_seq == 'multi')
      call baclose(lugb,igcret)
      call baclose(lugi,iicret)

    enddo ifhloop

------------------------------------------------------------------------
subroutine sort_storms_by_pressure
c
c   ABSTRACT: This subroutine  sorts storms by mslp.

  if (ifh > 1)
    call qsort


------------------------------------------------------------------------
subroutine get_next_ges

c   ABSTRACT: This subroutine calculates a guess position for the next
c         forecast time.

  call get_ij_bounds

  radmaxloop: for icut 1,icutmax .and. in_grid == 'n'
    levelloop: do n=1,nlevg
      if (readflag(ix1) .and. readflag(ix2)) ! the different level of wind

        call barnes
        
    enddo levelloop
  enddo radmaxloop

  call calcdist
  
------------------------------------------------------------------------
subroutine fixcenter

c   ABSTRACT: This subroutine loops through the different parameters
c         for the input storm number (ist) and calculates the 
c         center position of the storm by taking an average of
c         the center positions obtained for those parameters.



  do ip=1,maxtp
    if (calcparm(ip,ist)) ! test if parm definied (why ?)
      call calcdist

  if (iclose > 0) ! if there estimate of centers
    do ip=1,maxtp
      if (calcparm(ip,ist))
      call calcdist

    call avgcalc
    if (previous is success)
      call stdevcalc
      


    if (kprm > 0) ! what is kprm ??
      call wtavrg_lon
      call wtavrg


  if (itot4next > 0 .and. ifret /= 95) ! ??
    call stdevcalc
    
   
-----------------------------------------------------------------------
subroutine find_maxmin

c   This routine  finds the location (clon,clat) of and value of the
c   the max or min of fxy in the vicinity of slon,slat.


  if (cparm == 'vmag') then
  else
    call get_ij_bounds

  call date_and_time
  
  jloop: do j=-npts,npts,bskip1
  iloop: do i=-npts,npts,bskip1
    call check_valid_point
    call calcdist
    call barnes

  call get_ij_bounds 

  do k=1,nhalf !?
    call date_and_time
    jloop2: do j=-npts,npts,iskip
    iloop2: do i=-npts,npts,iskip
      call check_valid_point
      call barnes


------------------------------------------------------------------------
subroutine barnes

c   ABSTRACT: This routine performs a single-pass barnes anaylsis
c   of fxy at the point (flon,flat).

  do jix=jjbeg,jjend,bskip
  do iix=iibeg,iiend,bskip
    call calcdist

-----------------------------------------------------------------------
subroutine get_ij_bounds

c   ABSTRACT: This subroutine figures out, based on ri, dx and dy and
c   the guess latitude and longitude positions, the farthest reaching
c   grid points that are searchable by an analysis subroutine.


------------------------------------------------------------------------
subroutine check_valid_point

c   ABSTRACT: This subroutine checks to see if the input lat/lon
c   point is associated with four surrounding (i,j) locations that
c   have valid data. 

  call fix_latlon_to_ij 


------------------------------------------------------------------------
subroutine fix_latlon_to_ij

c   ABSTRACT: This subroutine takes an input lat/lon position and
c   assigns it to a nearby (i,j) gridpoint.



------------------------------------------------------------------------
subroutine first_ges_center

c   ABSTRACT: This subroutine scans an array and picks out areas of 
c   max or min, then loads those center positions into the first-
c   guess lat & lon arrays to be used by subroutine  tracker for 
c   locating the very specific low center positions.

  call find_all_maxmins
  
  

------------------------------------------------------------------------
subroutine find_all_maxmins

c   ABSTRACT: This subroutine will search an area delineated by  
c   input i and j indeces in order to find all local maxes or mins 
c   in that area.


  if (trkrinfo%choose_t2 == 'y') then
  else

    call avgcalc
    call stdevcalc
    
  search_loop: do while (still_finding_valid_maxmins)
    call get_ijplus1_check_wrap
    if (rough_gradient_check_okay) then

      call check_closed_contour
      
------------------------------------------------------------------------
subroutine check_closed_contour

c   ABSTRACT: This subroutine checks a field of data around an input
c   (ix,jx) data point to see if a closed contour exists around 
c   that data point.


  successive_contours_loop: do while (num_found_contours < num_requested_contours)
    single_contour_scan_loop: do while (still_scanning)
      multiple_ring_loop: do mr = 1,ringct
        call get_ijplus1_check_wrap

      individual_ring_loop: do ir = 1,9
        if (found_a_point_in_our_contour == 'y')
          if (get_last_isobar_flag == 'y')
            call calcdist

------------------------------------------------------------------------
SUBROUTINE qsort(x,ind,n)
c
c   Code converted using TO_F90 by Alan Miller
c   Date: 2002-12-18  Time: 11:55:47

   
   
