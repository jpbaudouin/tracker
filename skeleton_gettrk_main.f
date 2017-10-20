program trakmain
 !Main Program

    call date_and_time
    call w3tagb
    
    call read_nlists
    call read_fhours
    call read_tcv_card
    call read_gen_vitals

    if (inp%file_seq == 'onebig')
        if (inp%filetype == 'grib')
            call open_grib_files
        else
            call open_netcdf_files

    call tracker


    if (inp%filetype == 'grib')
        call baclose
        
    call w3tage
    
    
------------------------------------------------------------------------
subroutine tracker
    !Main loops
  
    if (inp%filetype == 'netcdf')
        call getgridinfo_netcdf
        
    ifhloop:  ifh = 1,ifhmax ! timesteps

        call date_and_time
        
        if (inp%file_seq == 'multi')
            call get_grib_file_name

            if (use_waitfor)
                call waitfor

            call open_grib_files
      
        if (inp%filetype == 'grib')
            call getgridinfo


        if (inp%filetype == 'grib') then
            call getdata
        else
          call getdata_netcdf


        !Computing Vorticity if needed
            call subtract_cor
            call rvcal
        !   


        if (trkrinfo%type == 'midlat' .or. trkrinfo%type == 'tcgen')
          call sort_storms_by_pressure
     

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
c          User wants us to run a command per forecast time

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
subroutine argreplace
      ! This subroutine is used to generate the pre-forecast-command


------------------------------------------------------------------------
subroutine open_grib_files

C     ABSTRACT: This subroutine must be called before any attempt is
C     made to read from the input GRIB files.  The GRIB and index files
C     are  opened with a call to baopenr.

    ! some conditions
        call baopenr
        call baopenw

------------------------------------------------------------------------
subroutine open_netcdf_files
c     Same as the previous one but for netcdf


------------------------------------------------------------------------
subroutine is_it_a_storm

c     ABSTRACT: This subroutine is called after the center of the storm
c     has been fixed.  Its purpose is to determine whether or not 
c     the center that was found is actually a storm, and not just some
c     passing trough (this has happened in the case of decaying or weak
c     storms).  It's called twice -- once to check for a minimum MSLP
c     gradient, and once to check for a circulation at 850 mb.  The 
c     subroutine input parameter "cparm" determines which parameter to
c     check for.


    call get_ij_bounds


    jloop: do jix = jbeg,jend,bskip
    iloop: do iix = ibeg,iend,bskip
  
        call calcdist

        if (defined_pt(i,j))

            if (cparm == 'v850')
                call getvrvt
              
    enddo iloop
    enddo jloop

------------------------------------------------------------------------
subroutine get_phase

c     ABSTRACT: This subroutine is a driver subroutine for
c     determining the structure or phase of a cyclone.  Initially, we
c     will just have it use the Hart cyclone phase space (CPS) scheme.


    if (phasescheme == 'cps' .or. phasescheme == 'both')
        if (ifh > 1)
            if (fixlon/lat initialised)
                call get_cps_paramb
                call get_cps_vth

    if (phasescheme == 'vtt' .or. phasescheme == 'both')
        call get_vtt_phase


------------------------------------------------------------------------
subroutine get_cps_paramb

c     ABSTRACT: This subroutine is part of the algorithm for determining
c     the structure, or phase, of a cyclone.  For Hart's cyclone phase
c     space, this subroutine determines "Parameter B", which determines
c     the degree of thermal symmetry between the "left" and "right" 
c     hemispheres of a storm, in the layer between 900 and 600 mb.


    call calcdist

    call get_ij_bounds


    jloop: do j=jbeg,jend,bskip
    iloop: do i=ibeg,iend,bskip

        call calcdist

    enddo iloop
    enddo jloop

------------------------------------------------------------------------
subroutine get_cps_vth
c     ABSTRACT: This subroutine is part of the algorithm for determining
c     the structure, or phase, of a cyclone.  For Hart's cyclone phase
c     space, this subroutine determines the thermal wind profile for 
c     either the lower troposphere (i.e., between 600 and 900 mb) or the
c     upper troposphere (i.e., between 300 and 600 mb)

    call get_ij_bounds
    

    levloop: do k = kbeg,kend
    jloop: do j=jbeg,jend,bskip
    iloop: do i=ibeg,iend,bskip

        call calcdist

    enddo iloop
    enddo jloop
    enddo levloop

    call calccorr

------------------------------------------------------------------------
subroutine calccorr

c     This subroutine is the main driver for a series of
c     other subroutines below this that will calculate the
c     correlation between two input arrays, xdat and ydat.

    call getmean

    call getdiff

    call getslope

    call getyestim
    
    call getresid

    call getcorr

------------------------------------------------------------------------
subroutine getmean
c     This subroutine is part of the correlation calculation,
c     and it simply returns the mean of the input array, xarr.

------------------------------------------------------------------------
subroutine getdiff
c     This subroutine is part of the correlation calculation,
c     and it returns in the array zdiff the difference values
c     between each member of the input array xarr and the
c     mean value, zmean.
     
------------------------------------------------------------------------
subroutine getslope
c     This subroutine is part of the correlation calculation,
c     and it returns the slope of the regression line.

------------------------------------------------------------------------
subroutine getyestim
c     This subroutine is part of the correlation calculation,
c     and it calculates all the predicted y-values using the
c     regression equation that has been calculated.

------------------------------------------------------------------------
subroutine getresid
c     This subroutine is part of the correlation calculation,
c     and it calculates all the residual values between the
c     input y data points and the y-estim predicted y values.

------------------------------------------------------------------------
subroutine getcorr
c     This subroutine is part of the correlation calculation,
c     and it does the actual correlation calculation.

-----------------------------------------------------------------------
subroutine get_vtt_phase
c     ABSTRACT: This subroutine is part of the algorithm for determining
c     the structure, or phase, of a cyclone.  Here, we are only looking
c     at the mid-to-upper tropospheric warm anomaly at the center of
c     the storm.


    call find_maxmin('tmp')

    if (previous call return initialised values)
        call fix_latlon_to_ij

    call check_closed_contour

------------------------------------------------------------------------
subroutine get_sfc_center
c     ABSTRACT: This subroutine computes a modified lat/lon fix position
c     to use as the input center position for the subroutines that
c     follow which calculate surface-wind related values. 

------------------------------------------------------------------------
subroutine get_wind_structure
c     ABSTRACT: This subroutine is a driver subroutine for
c     determining the structure of the low level winds of a cyclone

    if (ifh /= 1)
        call calcdist

    do iquad = 1,4
    do idist = 1,numdist

        call distbear
        call bilin_int_uneven
        if (previous is success)
            call getvrvt

    enddo
    enddo


    do iquad = 1,4
    do idist = 1,numdist

        call distbear
        call bilin_int_uneven
        if (previous is success)
            call getvrvt
            
------------------------------------------------------------------------
subroutine get_fract_wind_cov
c     ABSTRACT: This subroutine determines the fractional areal coverage
c     of winds exceeding various thresholds within specified arcs
c     (e.g., 200 km, 400 km, etc) in each quadrant of a storm.

    call get_ij_bounds


    jloop: do j = jbeg,jend
    iloop: do i = ibeg,iend

        call calcdist

            if (some conditions not realised) 
                njloop: do nj= ngridint,-ngridint,-1
                niloop: do ni= -ngridint,ngridint

                    call calcdist

            else

                call calcdist

------------------------------------------------------------------------
      subroutine get_ike_stats (imax,jmax,inp,dx,dy,ist,ifh
     &               ,fixlon,fixlat,xsfclon,xsfclat,valid_pt,calcparm
     &               ,ike,sdp,wdp,maxstorm,trkrinfo,igisret)
c
c     ABSTRACT: This subroutine computes the Integrated Kinetic Energy
c     (IKE) and Storm Surge Damage Potential (SDP) values, based on
c     Powell (BAMS, 2007).

    call get_ij_bounds
    
    do j = jbeg,jend
    do i = ibeg,iend

        call calcdist

    enddo
    enddo

------------------------------------------------------------------------
subroutine distbear
c
c     ABSTRACT: Given an origin at latitude, longitude=xlato,xlono,
c     this subroutine will locate a target point at a distance dist in
c     km or nautical miles (depends on what you use for "rad_earth..."
c     below), at bearing bear (degrees clockwise from north).
c     Returns latitude xlatt and longitude xlont of target point.
c

------------------------------------------------------------------------
subroutine bilin_int_uneven
c
c     ABSTRACT: This subroutine performs a bilinear interpolation to get
c     a data value at a given lat/lon that may be anywhere within a box
c     defined by the four surrouding grid points.

------------------------------------------------------------------------
subroutine sort_storms_by_pressure
c
c     ABSTRACT: This subroutine  sorts storms by mslp.

    if (ifh > 1)
        call qsort

------------------------------------------------------------------------
subroutine getvrvt
c
c     ABSTRACT: This subroutine takes as input a u-wind and v-wind value
c     at an input (xlon,xlat) location and returns the tangential and
c     radial wind components relative to the input center lat/lon 
c     position (centlon,centlat).


    call calcdist

    if (xlondiff == 0 .and. xlatdiff > 0)
    else if (xlondiff > 0 .and. xlatdiff == 0)
    else
        call calcdist


------------------------------------------------------------------------
subroutine output_atcfunix

c     ABSTRACT: This subroutine  outputs a 1-line message for a given 
c     storm at an input forecast hour in the new ATCF UNIX format.  



------------------------------------------------------------------------
subroutine output_hfip

c     ABSTRACT: This subroutine  outputs a 1-line message for a given
c     storm at an input forecast hour in a modified ATCF UNIX format.


------------------------------------------------------------------------
subroutine output_fract_wind

c     ABSTRACT: This subroutine  outputs a 1-line message for a given
c     storm at an input forecast hour.  This message contains the
c     values for the fractional areal coverage of various wind
c     thresholds.  In addition, this subroutine also writes out
c     records to a file containing data on the PDF of wind magnitudes
c     within r=350 km.
c

------------------------------------------------------------------------
subroutine output_wind_structure

c     ABSTRACT: This subroutine  outputs a 1-line message for a given
c     storm at an input forecast hour.  This message contains the
c     values of the winds at specified distances along 45-degree
c     radials in each storm quadrant.

------------------------------------------------------------------------
subroutine output_ike

c     ABSTRACT: This subroutine  outputs a 1-line message for a given
c     storm at an input forecast hour.  This message contains the values
c     for the Integrated Kinetic Energy (IKE) and Storm Surge Damage
c     Potential (SDP), based on Powell (BAMS, 2007).


------------------------------------------------------------------------
subroutine output_phase 

c     ABSTRACT: This subroutine  outputs a 1-line message for a given
c     storm at an input forecast hour.  This message contains the values
c     for the three parameters that comprise Bob Hart's cyclone phase
c     space (CPS).


------------------------------------------------------------------------
subroutine output_atcf_gen

c     ABSTRACT: This subroutine  outputs a 1-line message for a given 
c     storm at an input forecast hour in a modified atcfunix format.  


------------------------------------------------------------------------
subroutine output_atcf_sink

c     ABSTRACT: This subroutine  outputs a 1-line message for a given 
c     storm at an input forecast hour in a modified atcfunix format.  


------------------------------------------------------------------------
subroutine output_tcvitals

c     ABSTRACT: This subroutine  outputs a tcvitals record.


----------------------------------------------------------------------
subroutine output_gen_vitals

c     ABSTRACT: This subroutine  outputs a modified vitals record.  


------------------------------------------------------------------------
subroutine get_next_ges

c     ABSTRACT: This subroutine calculates a guess position for the next
c               forecast time.

    call get_ij_bounds

    radmaxloop: for icut 1,icutmax .and. in_grid == 'n'
        levelloop: do n=1,nlevg
            if (readflag(ix1) .and. readflag(ix2)) ! the different level of wind

                call barnes
                
        enddo levelloop
    enddo radmaxloop

    call calcdist
    

------------------------------------------------------------------------
subroutine getradii
c
c     ABSTRACT: This subroutine looks through the wind data near an
c     input storm center (fixlon,fixlat) and gets the radii of various
c     surface winds in each of the 4 storm quadrants (NE,NW,SE,SW).  

    jloop: do j=jbeg,jend
    iloop: do i=ibeg,iend

        call calcdist


    call qsort



------------------------------------------------------------------------
subroutine get_max_wind
      
c     ABSTRACT: This subroutine looks for the maximum near-surface wind
c     near the storm center.


    do j=jbeg,jend
    do i=ibeg,iend

        call calcdist
        
------------------------------------------------------------------------
subroutine fixcenter

c     ABSTRACT: This subroutine loops through the different parameters
c               for the input storm number (ist) and calculates the 
c               center position of the storm by taking an average of
c               the center positions obtained for those parameters.



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
        
        
------------------------------------------------------------------------
subroutine avgcalc

c     ABSTRACT: This subroutine just calculates a straight average of
c     the parameters in the input array (xdat). 


------------------------------------------------------------------------
subroutine wtavrg
c
c     ABSTRACT: This subroutine calculates a weighted average of the 
c     parameters in the input array (xdat) using the input weights
c     in the input array (wt). 


-----------------------------------------------------------------------
subroutine wtavrg_lon

c     ABSTRACT: This subroutine calculates a weighted average of the
c     parameters in the input array (xlon) using the input weights
c     in the input array (wt).


-----------------------------------------------------------------------
subroutine stdevcalc

! ??


------------------------------------------------------------------------
subroutine get_uv_center

c     ABSTRACT: This subroutine calculates the center fix position for
c     the minimum in the wind speed near the storm center.


    call get_ij_bounds
    
    


    do intnum = 1,numinterp ! ??

        call bilin_int_even
        call lin_int_lon
        call lin_int
 
 
    call calc_vmag

    call find_maxmin


------------------------------------------------------------------------
subroutine get_uv_guess

c     ABSTRACT: The purpose of this subroutine is to get a modified 
c               first guess lat/lon position before searching for the 
c               minimum in the wind field.


    do ip = 1,maxtp
        if ((ip > 2 .and. ip < 7) .or. ip == 10) then
        else
          if (calcparm(ip,ist)) then
            call calcdist

------------------------------------------------------------------------
subroutine calc_vmag

c     ABSTRACT: This subroutine calculates the magnitude of the wind
c     speed for an array of points, given real u and real v arrays.

------------------------------------------------------------------------
subroutine bilin_int_even

c     ABSTRACT: This subroutine does a bilinear interpolation on a 
c     grid of evenly spaced data.


------------------------------------------------------------------------
subroutine lin_int
c
c     ABSTRACT: This subroutine linearly interpolates evenly spaced
c               data from one grid to another.
c 

------------------------------------------------------------------------
subroutine lin_int_lon

c     ABSTRACT: This subroutine linearly interpolates evenly spaced
c               data from one grid to another.


------------------------------------------------------------------------
subroutine get_zeta_values
c
c     ABSTRACT: This subroutine finds the maximum and mean zeta values
c     at 850 & 700 mb, near a storm center. 

    call get_ij_bounds
        report_zeta_loop: do n=1,2
            if (zeta(ilonfix,jlatfix,n) > -9990.0)
                call find_maxmin
                
            call fix_latlon_to_ij
            
-----------------------------------------------------------------------
subroutine find_maxmin

c     This routine  finds the location (clon,clat) of and value of the
c     the max or min of fxy in the vicinity of slon,slat.


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

c     ABSTRACT: This routine performs a single-pass barnes anaylsis
c     of fxy at the point (flon,flat).

    do jix=jjbeg,jjend,bskip
    do iix=iibeg,iiend,bskip
        call calcdist

-----------------------------------------------------------------------
subroutine get_ij_bounds

c     ABSTRACT: This subroutine figures out, based on ri, dx and dy and
c     the guess latitude and longitude positions, the farthest reaching
c     grid points that are searchable by an analysis subroutine.


------------------------------------------------------------------------
subroutine check_bounds

------------------------------------------------------------------------
subroutine calcdist
c
c     ABSTRACT: This subroutine computes the distance between two 
c               lat/lon points by using spherical coordinates to 
c               calculate the great circle distance between the points.


------------------------------------------------------------------------
subroutine subtract_cor
c
c     ABSTRACT: This subroutine  subtracts out the coriolis parameter
c     from the vorticity values.

------------------------------------------------------------------------
subroutine get_grib_file_name

c     ABSTRACT: This subroutine uses various input regarding the model
c     and forecast hour and generates the name of the input grib file
c     for this particular forecast hour.


------------------------------------------------------------------------
subroutine getdata
      
c     ABSTRACT: This subroutine reads the input GRIB file for the
c     tracked parameters. 

    do ip = 1,nparms
        call getgb
        if (previous is success)
            call bitmapchk
            if (lbrdflag .eq. 'n')
                call conv1d2d_logic

            select case (chparm(ip))
                call conv1d2d_real
        else


    if (phaseflag == 'y') then
        if (phasescheme == 'cps' .or. phasescheme == 'both') then
            do ip = 1,nlevs_cps
                call getgb

            if (previous is success)
                call bitmapchk
                call conv1d2d_real

------------------------------------------------------------------------
subroutine bitmapchk
c
c     This subroutine checks the bitmap for non-existent data values.


-----------------------------------------------------------------------
subroutine getdata_netcdf

c     ABSTRACT: This subroutine reads the input netcdf file for the
c     tracked parameters.
c

    call find_var ! numerous call for different variables with different context
    

-----------------------------------------------------------------------
subroutine find_var

c     ABSTRACT: This subrtoutine reads the input of a netcdf file,
c     tests if a variable have one of the given names, and stock
c     its values.


------------------------------------------------------------------------
subroutine conv1d2d_logic

c     ABSTRACT: This subroutine converts a 1-dimensional input 
c     array of logical data (lb1d) into a 2-dimensional output
c     array (dimension imax,jmax) of logical data (lb2d).


------------------------------------------------------------------------
subroutine conv1d2d_real
c
c     ABSTRACT: This subroutine converts a 1-dimensional input 
c     array of real data (dat1d) into a 2-dimensional output
c     array (dimension imax,jmax) of real data (dat2d).
c

------------------------------------------------------------------------
subroutine read_nlists

c     ABSTRACT: This subroutine simply reads in the namelists that are
c     created in the shell script. 


------------------------------------------------------------------------
subroutine read_fhours

c     ABSTRACT: This subroutine reads in a text file that contains the
c     forecast times that will be read in.

---------------------------------------------------------------------
subroutine read_tcv_card

c     ABSTRACT: This subroutine reads in the updated TC Vitals file
c               for the current time and prints out those cards (storms)
c               that have been selected to be processed. 


------------------------------------------------------------------------
subroutine read_gen_vitals

c     ABSTRACT: This subroutine reads in a modified TC Vitals file
c     for the current time and prints out those cards (storms) that
c     have been selected to be processed.


------------------------------------------------------------------------
subroutine getgridinfo

c     ABSTRACT: The purpose of this subroutine is just to get the max
c     values of i and j and the dx and dy grid spacing intervals for the
c     grid to be used in the rest of the program.

    call getgb


------------------------------------------------------------------------
subroutine getgridinfo_netcdf


c     ABSTRACT: The purpose of this subroutine is just to get the max
c     values of i and j and the dx and dy grid spacing intervals for the
c     grid to be used in the rest of the program. 


------------------------------------------------------------------------
subroutine check_valid_point

c     ABSTRACT: This subroutine checks to see if the input lat/lon
c     point is associated with four surrounding (i,j) locations that
c     have valid data. 

    call fix_latlon_to_ij 


------------------------------------------------------------------------
subroutine fix_latlon_to_ij

c     ABSTRACT: This subroutine takes an input lat/lon position and
c     assigns it to a nearby (i,j) gridpoint.



------------------------------------------------------------------------
subroutine rvcal

c     ABSTRACT: This routine calculates the relative vorticity (zeta)
c     from u,v on an evenly-spaced lat/lon grid.


------------------------------------------------------------------------
subroutine first_ges_center

c     ABSTRACT: This subroutine scans an array and picks out areas of 
c     max or min, then loads those center positions into the first-
c     guess lat & lon arrays to be used by subroutine  tracker for 
c     locating the very specific low center positions.

    call find_all_maxmins
    
  

------------------------------------------------------------------------
subroutine find_all_maxmins

c     ABSTRACT: This subroutine will search an area delineated by  
c     input i and j indeces in order to find all local maxes or mins 
c     in that area.


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

c     ABSTRACT: This subroutine checks a field of data around an input
c     (ix,jx) data point to see if a closed contour exists around 
c     that data point.


    successive_contours_loop: do while (num_found_contours < num_requested_contours)
        single_contour_scan_loop: do while (still_scanning)
            multiple_ring_loop: do mr = 1,ringct
                call get_ijplus1_check_wrap

            individual_ring_loop: do ir = 1,9
                if (found_a_point_in_our_contour == 'y')
                    if (get_last_isobar_flag == 'y')
                        call calcdist


------------------------------------------------------------------------
subroutine get_ijplus1_check_wrap

c     ABSTRACT: This subroutine takes an (i,j) position input and 
c     returns the four neighboring (i,j) points to the east, south, 
c     west and north. 

------------------------------------------------------------------------
SUBROUTINE qsort(x,ind,n)
c
c     Code converted using TO_F90 by Alan Miller
c     Date: 2002-12-18  Time: 11:55:47
