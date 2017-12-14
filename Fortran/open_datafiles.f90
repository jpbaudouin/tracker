!#######################################################################
!######### DEFINITION OF THE SUBROUTINES TO READ THE DATAFILES #########
!#######################################################################


!25-10-2017 Creation
!12-12-2017 1st Commit -> 1st version

! RQ for further developments:
! Is it possible that a the time dimension is not oredered in a netcdf?
! Is ts_num really useful??
! change the scope of module datafile_param according to what is actually
! used in the main program and in the other subroutines
! I need to keep the number of global variables to a minimum, the one very often used
! Check if data_var needs explicit shape in find_var
! Add comments in getgridinfo and getdata about what has been read


! CONTAINS :
!   The module "data_param" with all the variable realtive to the data
! 
!   The module "datafile" with the distributive subroutines :
!       - open_file
!       - getgridinfo
!       - getdata(ts)
! 
! 
! 
!   The module "netcdf_datafile" for the reading of the data from a netcdf file,
!   It contains the subroutines:
!       - open_netcdf_file
!       - getgridinfo_netcdf
!       - getdata_netcdf(ts)
!       - find_var(nvar, name_var, level, ts, readflag, data_var),
!           called by getdata_netcdf
! 
!   The module "grib_datafile" for the reading of the data from a netcdf file,
!   It contains the subroutines:
!       - open_grib_file
!       - getgridinfo_grib
!       - getdata_grib
!   !!!!!!!!!!!!!!!! --> Not yet implemented see module_grid_datafiles.f90.draft
! 
! 
! 
!   The module "time_step" including the subroutines:
!       - readts, calling either the two following :
!       - readfort15
!       - readts_netcdf


!#######################################################################
!############ Module containing all the relevant variables #############
!######## concerning the datafiles (timestep, grid and fields) #########
!#######################################################################
!-----------------------------------------------------------------------

module data_param

! imax        integer number of pts in i-direction on grid
! jmax        integer number of pts in j-direction on grid
! dx       grid spacing in i-direction
! dy       grid spacing in j-direction
! valid_pt Logical; bitmap indicating if valid data at that pt
  
  
  implicit none
  
  type readflags  
    logical detection
    logical, allocatable :: fixcenter(:)
    logical, allocatable :: advection(:)
    logical, allocatable :: intensity(:)
  end type readflags
    
  type (readflags) Rflag

  ! from module grid_bounds
  real, save ::  glatmin, glatmax, glonmin, glonmax  ! These 
                 ! define the boundaries of the input data grid
  real, save, allocatable ::  glat(:), glon(:)  ! Will be filled
                 ! with lat/lon values for each pt on input grid
                 
  integer, save :: imax, jmax
  real, save :: dx, dy
  
  real, save, allocatable :: data_detection(:,:)
  real, save, allocatable :: data_fixcenter(:,:,:)
  real, save, allocatable :: data_intensity(:,:,:)
  
  real, save, allocatable :: data_u(:,:,:)
  real, save, allocatable :: data_v(:,:,:)
  

     


end module data_param


 
  
  
!#################################################################################
!################# Module to read a datafile using netcdf format #################
!#################################################################################
!---------------------------------------------------------------------------------

module netcdf_datafile

  use namelist_trk
  use data_param
  use netcdf
  
  implicit none
  
  integer, save :: ncid
  logical, save :: need_to_flip_lats, need_to_flip_lons
  
  ! need_to_flip_lats logical flag read in from getgridinfo that
  !   indicates if data needs flipped north to south
  ! need_to_flip_lons logical flag read in from getgridinfo that
  !   indicates if data needs flipped east to west
  ! ncid is the id of the datafile

  contains
  
  !---------------------------------------------------------------------------------
  !--------------------- Subroutine opening a netcdf datafile ----------------------
  !---------------------------------------------------------------------------------  
  subroutine open_netcdf_file
  ! This subroutine reproduce the behaviour of open_grib_files
  ! but is not as important, it sets the id (ncid) of the opened file
  
    implicit none
    integer ret
  
    ret = nf90_open("fort.11",nf90_nowrite,ncid)
  
    if (ret /= nf90_noerr) then
      if ( verb .ge. 1 ) then
        print *,' '
        print *,'!!! ERROR 130 in subroutine open_netcdf_file : '
        print *,'!!! nf90_open return codes :', ret
      endif
      STOP 130
    endif
    
  end subroutine
  
  !---------------------------------------------------------------------------------
  !---------- Subroutine getting grid information from a netcdf datafile -----------
  !---------------------------------------------------------------------------------    
  subroutine getgridinfo_netcdf
    !ABSTRACT: The purpose of this subroutine is just to get the max
    ! values of i and j and the dx and dy grid spacing intervals for the
    ! grid to be used in the rest of the program.  So just read the 
    ! netcdf file to get the lon and lat data.  Also, get the info for 
    ! the data grid's boundaries.  This boundary information will be 
    ! used later in the tracking algorithm, and is accessed via Module 
    ! grid_bounds.
   
    implicit none
  
    integer  nDim,nVar,nAtt,unlimitDimId,formatNum
    integer ret1,ret2,ret3, erra_a1, erra_a2
    integer idlat, idlon, i
    character(len=nf90_max_name) namelat,namelon


    !-------------------------------------------------------------------
    ! We first find the dimensions for latitude (jmax) and 
    ! longitude (imax) using inquiry on the dimension of 
    ! the nectdf file
    !-------------------------------------------------------------------
  
    ret1 = nf90_inquire (ncid,nDim,nVar,nAtt,unlimitDimId,formatNum)
    if (ret1 /= NF90_NOERR) then
      if ( verb .ge. 1 ) then
        print *,' '
        print *,'!!! ERROR 131 in getgridinfo_netcdf reading file'
        print *,'!!! nf90_inquire ret1 = ', ret1
      endif
      STOP 131
    endif
  

    ! scanning latitude / longitude name
    namelat = "latitude"
    ret1 = nf90_inq_dimid(ncid, namelat, idlat)
    if (ret1 /= NF90_NOERR) then
      namelat = "lat"
      ret2 = nf90_inq_dimid(ncid, namelat, idlat)
      if (ret2 /= NF90_NOERR) then
        namelat = "y"
        ret3 = nf90_inq_dimid(ncid, namelat, idlat)
        if (ret3 /= NF90_NOERR) then
          if ( verb .ge. 1 ) then
            print *,' '
            print *,'!!! ERROR 132 in getgridinfo_netcdf when finding'
            print *,'!!! the latitude, call to nf90_inq_dimid :'
            print *,'!!! inquiry name : "latitude", ret1 =',ret1
            print *,'!!! inquiry name : "lat", ret2 =',ret2
            print *,'!!! inquiry name : "y", ret3 =',ret3
          endif
          STOP 132
        endif
      endif
    endif
    
    namelon = "longitude"
    ret1 = nf90_inq_dimid(ncid, namelon, idlon)
    if (ret1 /= NF90_NOERR) then
      namelon = "lon"
      ret2 = nf90_inq_dimid(ncid, namelon, idlon)
      if (ret2 /= NF90_NOERR) then
        namelon = "x"
        ret3 = nf90_inq_dimid(ncid, namelon, idlon)
        if (ret3 /= NF90_NOERR) then
          if ( verb .ge. 1 ) then
            print *,' '
            print *,'!!! ERROR 133 in getgridinfo_netcdf when finding'
            print *,'!!! the longitude, call to nf90_inq_dimid :'
            print *,'!!! inquiry name : "longitude", ret1 =',ret1
            print *,'!!! inquiry name : "lon", ret2 =',ret2
            print *,'!!! inquiry name : "x", ret3 =',ret3
          endif
          STOP 132
        endif
      endif
    endif
  
  
        
    ret1 = nf90_inquire_dimension (ncid, idlon, namelon, imax)
    ret2 = nf90_inquire_dimension (ncid, idlat, namelat, jmax)
    if (ret1 /= NF90_NOERR .or. ret2 /= NF90_NOERR) then
      if ( verb .ge. 1 ) then
        print *,' '
        print *,'!!! ERROR in 134 getgridinfo_netcdf finding dimensions'
        print *,'!!! inquiry longitude dimension, ret1 =',ret1
        print *,'!!! inquiry latitude  dimension, ret2 =',ret2
      endif
      STOP 134
    endif   
  
    if (allocated(glat)) deallocate(glat)
    if (allocated(glon)) deallocate(glon)
    allocate (glon(imax), stat = erra_a1)
    allocate (glat(jmax), stat = erra_a2)
    if (erra_a1 /= 0 .or. erra_a2 /= 0) then
      if ( verb .ge. 1 ) then
        print *,' '
        print *,'!!! ERROR 135 in getgridinfo_netcdf allocating '
        print *,'!!! glat : imax = ', imax, ", err = ", erra_a1
        print *,'!!! glon : jmax = ', jmax, ", err = ", erra_a2
      endif
      STOP 135
    endif
        
    if ( verb .ge. 3 ) then
      print *,' '
      print "(A, A15, A, I5)",'longitude name: ', trim(namelon), &
              ', longitude dimension: ', imax
      print "(A, A15, A, I5)",'latitude  name: ', trim(namelat), &
              ', latitude  dimension: ', jmax
    endif
  
  
  
    !-------------------------------------------------------------------
    ! Then, we fill the other grid parameters : glon and glat for the
    ! values, dx and dy for the distance between each grid point),
    ! using inquiry on the variables of the netcdf file  
    !-------------------------------------------------------------------

    ret1 = nf90_inq_varid(ncid,namelon,idlon)
    ret2 = nf90_inq_varid(ncid,namelat,idlat)
    if (ret1 /= NF90_NOERR .or. ret2 /= NF90_NOERR) then
      if ( verb .ge. 1 ) then
        print *,' '
        print *,'!!! ERROR 136 in getgridinfo_netcdf :'
        print *,'!!! dimensions were not found among variables'
        print *,'!!! inquiry longitude varid, ret1 =',ret1
        print *,'!!! inquiry latitude  varif, ret2 =',ret2
      endif
      STOP 136
    endif  
      
    ret1 = nf90_get_var(ncid, idlon, glon)
    ret2 = nf90_get_var(ncid, idlat, glat)
    if (ret1 /= NF90_NOERR .or. ret2 /= NF90_NOERR) then
      if ( verb .ge. 1 ) then
        print *,' '
        print *,'!!! ERROR 137 in getgridinfo_netcdf '
        print *,'!!! when filling dimensions values'
        print *,'!!! inquiry longitude values, ret1 =',ret1
        print *,'!!! inquiry latitude  values, ret2 =',ret2
      endif
      STOP 137
    endif
      
    if (minval(glat) < -90 .or. maxval(glat) > 90) then
      if ( verb .ge. 1 ) then
        print *,'!!! ERROR 138 in getgridinfo_netcdf '
        print *,'!!! on the boundary of the latitude grid :'
        print *,'!!! minimum latitude :',minval(glat)
        print *,'!!! maximum latitude :',maxval(glat)
      endif
    STOP 138
    endif
      
    if ( (minval(glon) < -180 .or. maxval(glon) > 540) .or. &
         (minval(glon) < 0   .and. maxval(glon) > 180)) then
      if ( verb .ge. 1 ) then
        print *,'!!! ERROR 139 in getgridinfo_netcdf '
        print *,'!!! on the boundary of the longitude grid :'
        print *,'!!! minimum longitude :',minval(glon)
        print *,'!!! maximum longitude :',maxval(glon)
      endif
    STOP 139
    endif
      
    if (minval(glon(1:(imax-1))-glon(2:imax)) /= &
        maxval(glon(1:(imax-1))-glon(2:imax))) then
      if (minval(glat(1:(imax-1))-glat(2:imax)) /= &
          maxval(glat(1:(imax-1))-glat(2:imax))) then
        if ( verb .ge. 1 ) then
          print *,' '
          print *,'!!! ERROR 140 in getgridinfo_netcdf when reading'
          print *,'!!! dimensions values : the lat/lon steps between'
          print *,'!!! each point is not constant.'
          print *,'!!! Gaussian grid unsupported'
        endif
        STOP 140
      endif 
    endif
      
    dx = glon(2)-glon(1)
    dy = glat(2)-glat(1)
      
    need_to_flip_lons = .false.
    if (dx < 0) then
      dx = -dx
      need_to_flip_lons = .true.
      glon = glon(imax:1:-1)
    endif

    need_to_flip_lats = .false.
    if (dy < 0) then
      dy = -dy
    else
      need_to_flip_lats = .true.
      glat = glat(jmax:1:-1)
    endif 
   
    
    !-------------------------------------------------------------------
    ! Finally, we fill the value for the grid boundaries and check
    ! to see if the requested boundary limits that the user input
    ! the user input are contained within this grid
    !-------------------------------------------------------------------

    glatmin = minval(glat)
    glatmax = maxval(glat)
      
    glonmin = minval(glon)
    glonmax = maxval(glon)
    
    if ( verb .ge. 3 ) then
      print *,' '
      print *,'Data Grid Lat/Lon boundaries follow:'
      print "(A, f8.3, A, f8.3)", 'Min Lat: ', glatmin, &
                                 ' Min Lon: ', glonmin
      print "(A, f8.3, A, f8.3)", 'Max Lat: ', glatmax, &
                                 ' Max Lon: ', glonmax
    endif
    
    if (fileinfo%grid_type == 'regional') then
      if (trkrinfo%eastbd > glonmax) then
        trkrinfo%eastbd = glonmax - 5.0
        go to 141
      endif   

      if (trkrinfo%westbd < glonmin) then
        trkrinfo%westbd = glonmin + 5.0
        go to 141
      endif 
          
      if (trkrinfo%northbd > glatmax) then
        trkrinfo%northbd = glatmax - 5.0          
        go to 141
      endif 
       
      if (trkrinfo%southbd < glatmin) then
        trkrinfo%southbd = glatmin + 5.0            
        go to 141
      endif 
    
      return
    
      141 continue 
      if ( verb .ge. 3 ) then
        print *, " "
        print *, "WARNING the user requested bourndaries are beyond"
        print *, "the boundaries of the data as defenied in the datafile."
        print *, "The tracker boundaries are therefore modified"
      endif
    endif
  
  end subroutine  
  


  !---------------------------------------------------------------------------------
  !------------ Subroutine getting fields data from a netcdf datafile --------------
  !---------------------------------------------------------------------------------    
  subroutine getdata_netcdf(ts)
    ! ABSTRACT: This subroutine reads the input netcdf file for the
    ! tracked parameters.
    
    ! ts was previously ifh, it is now the same value as found in the datafile (ts_tim(its))
    
    ! detection (1)
    ! fixcenter (")
    ! u / v (advection)
    ! intensity (")

    integer, intent(in) :: ts
    
    logical Rflag_adv1(numfield%advection)
    logical Rflag_adv2(numfield%advection)
    logical Rflag_detection(1)
    integer err_a1, err_a2
    real data_detection2(imax, jmax,1)
 
    
    ! loading the field for detection of new lows
    if (trkrinfo%genesis_flag) then 
      call find_var(1, [fname%detection], [fname%detection_lev], ts, &
                    Rflag_detection, data_detection2)
      Rflag%detection = Rflag_detection(1)
      data_detection  = data_detection2(:,:,1)
      if (.not. Rflag%detection) then
        if (verb .ge. 1) then
          print *, ""
          print *, "!!! ERROR 151: in getdata_netcdf"
          print *, "!!! The field ", fname%detection, "wasn't found"
          print *, "!!! in the netcdf file for the timestep ", ts
          print *, "!!! while trkrinfo%genesis_flag is true"
        endif 
        STOP 151
      endif                    
    else
      Rflag%detection = .false.
    endif


    ! loading the fields for fixing the center    
    call find_var(numfield%fixcenter, fname%fixcenter, &
                  fname%fixcenter_lev, ts, &
                  Rflag%fixcenter, data_fixcenter)  

    if (.not. any(Rflag%fixcenter)) then
      if (verb .ge. 1) then
        print *, ""
        print *, "WARNING: in getdata_netcdf, not all the fixcenter"
        print *, "fields where founds for the timestep", ts
        print *, "fname%fixcenter names: ", fname%fixcenter
        print *, "fname%fixcenter read: ", Rflag%fixcenter
        print *, "The processing will continue with the available data"
      endif 
    endif    


    ! loading the fields for the center advection
    call find_var(numfield%advection, fname%adv_1stdim, &
                  fname%adv_lev, ts, Rflag_adv1, data_u)
                  
    call find_var(numfield%advection, fname%adv_2nddim, &
                  fname%adv_lev, ts, Rflag_adv2, data_v)  
    
    Rflag%advection = Rflag_adv1 .and. Rflag_adv2
    if (.not. any(Rflag%advection)) then
      if (any(.not. Rflag%advection)) then
        if (verb .ge. 1) then
          print *, ""
          print *, "!!! ERROR 153 : in getdata_netcdf"
          print *, "!!! No advection fields where founds in the netcdf"
          print *, "!!! file for the timestep", ts
          print *, "!!! fname%adv_1stdim names: ", fname%adv_1stdim
          print *, "!!! fname%adv_2nddim names: ", fname%adv_2nddim
        endif
        STOP 153
      else
        if (verb .ge. 1) then
          print *, ""
          print *, "WARNING: in getdata_netcdf, not all the advection"
          print *, "fields where founds for the timestep", ts
          print *, "fname%adv_1stdim names: ", fname%adv_1stdim
          print *, "fname%adv_2nddim names: ", fname%adv_2nddim
          print *, "fname%advection read: ", Rflag_adv1, Rflag_adv2
          print *, "The tracking will continue with the available data"
        endif
      endif
    
    endif      

    
    ! loading the fields for the center advection
    call find_var(numfield%intensity, fname%intensity,  &
                  fname%intensity_lev, ts, &
                  Rflag%intensity, data_intensity)      

    if (.not. any(Rflag%intensity)) then
      if (verb .ge. 1) then
        print *, ""
        print *, "WARNING: in getdata_netcdf, not all the intensity"
        print *, "fields where founds for the timestep", ts
        print *, "fname%intensity names:", fname%intensity
        print *, "fname%intensity read", Rflag%intensity
        print *, "The processing will continue with the available data"
      endif 
    endif 
    
    
  end subroutine

  
  !---------------------------------------------------------------------------------
  !------------ Subroutine getting fields data from a netcdf datafile --------------
  !---------------------------------------------------------------------------------    
  subroutine find_var(nvar, name_var, level, ts, readflag, data_var)
  
  ! ABSTRACT: This subrtoutine reads the input of a netcdf file,
  ! tests if a variable have one of the given names, and stock
  ! its values. The variable needs to have 4 or 3 dimensions.
  
  ! INPUT:
  ! nvar      number of variables (dimension of name_var, level,
  !                                readflag, and first of data_var)
  ! name_var  name  of the variables to be searched
  ! level     level of the variables
  ! ts        timestep value as find in the datafile
  
  ! OUTPUT
  ! readflag  logical array, indicates if a parm was read in
  ! data_var  array which stock the value of the variable
  
  
  ! OTHER:
  ! id_dim    vector of dimension 4 giving the place of each 
  !           dimension for the array read in the netcdf files,
  !           in the order : longitude, latitude, time, level  

 
    implicit none
  
    integer, intent(in) :: nvar, ts
    integer, intent(in) :: level(:) !nvar
    character(len=16), intent(in) :: name_var(:) !nvar
    real,    intent(out) :: data_var(:,:,:) ! imax, jmax, nvar
    logical, intent(out) :: readflag(:)  !nvar
  
    character(len=NF90_MAX_NAME) namedim
    real, allocatable :: values (:,:)
    integer varid,levelid,timeid,ndims,length
    real missingValue
    integer dimids(NF90_MAX_VAR_DIMS)
    integer,allocatable :: valueslvl(:),valuestime(:)
    integer iddim(4),starts(4),counts(4)
    integer n,id,t,l
    integer ret1,ret2,ret3,ret4,ret5,ret6



    iddim = (/0,0,0,0/) ! (lon,lat,time,level) 
  
    if ( verb .ge. 3 ) then
      print *,' '
      print *,'Beginning of the subroutine find_var'
    endif 
    
    
    nameloop : do n = 1,nvar
    
      ret1 = nf90_inq_varid(ncid,name_var(n),varid)
      if (ret1 == NF90_NOERR) then      
        ret2 = nf90_inquire_variable(ncid, varid, ndims = ndims, &
                                                 dimids = dimids)
      endif
      
      if (ret1 /= NF90_NOERR .or. ret2 /= NF90_NOERR) then
        if (verb .ge. 1) then
          print *, ""
          print *, "!!! ERROR 161 in find_var"
          print *, "!!! the variable ", name_var(n)
          print *, "!!! isn't found or can't be inquired and is skipped"
        endif
        go to 160 ! skip name
      endif


    !-----------------------
    ! Testing the dimension:
    !  We suppose that the level doesn't have the same length
    !  as the latitude or the longitude, and that the time 
    !  dimension is called 'time'.
    !--------------------------------------------------------   
      if (ndims == 4 .or. ndims == 3) then
        do id = 1,ndims
          ret3 = nf90_inquire_dimension(ncid,dimids(id),namedim,length)

          if (ret3 == NF90_NOERR) then
            if (namedim == 'time') then
              if (iddim(3) == 0) then
                iddim(3) = id
                go to 171 ! continue the analysis
              endif               
            else if (length == imax) then
              if (iddim(1) == 0) then
                iddim(1) = id
                go to 171 ! continue the analysis
              endif
            else if (length == jmax) then
              if (iddim(2) == 0) then
                iddim(2) = id
                go to 171 ! continue the analysis
              endif
            else 
              if (iddim(4) == 0 .and. ndims == 4) then
                iddim(4) = id
                go to 171 ! continue the analysis
              endif
            endif
          endif

          if ( verb .ge. 1 ) then
            print *,' '
            print *, '!!! Error 162 in subroutine find_var'
            print *, '!!! when analysing the dimensions'
            print *, '!!! of the variable ', name_var(n)
            print *, '!!! The variable is skipped'
          endif
          go to 160 ! skip name
                
          171 continue !from a good analysis of a dimension
                
        enddo

      else
        if (verb .ge. 1) then
          print *, ""
          print *, "!!! ERROR 163 in find_var"
          print *, "!!! the variable detected ", name_var(n)
          print *, "!!! has a wrong number of dimensions: ", ndims
          print *, "!!! 3 or 4 are expected (lon, lat, 'time', -level-)"
          print *, "!!! the variable is skipped"
        endif
        go to 160 ! skip name
      endif
         
            
    ! inquiring the values for the level dimension if require
    !--------------------------------------------------------
      if (ndims == 4) then
        ret3 = nf90_inquire_dimension(ncid,dimids(iddim(4)), &
                                      namedim,length)
        ret4 = nf90_inq_varid(ncid,namedim,levelid)
        if (ret3 == NF90_NOERR .and. ret4 == NF90_NOERR) then
          if (allocated(valueslvl)) deallocate(valueslvl)
          allocate (valueslvl(length),stat=ret5)
          if (ret5 == 0) then
            ret6 = nf90_get_var(ncid,levelid,valueslvl)
            if (ret6 == NF90_NOERR) then
              do l = 1,length
                if(valueslvl(l) == level(n)) then
                  go to 172 ! good analysis, l is the level required
                endif
              enddo
            endif
          endif
        endif
      endif
      
      if ( verb .ge. 1 ) then
        print *,' '
        print *, '!!! Error 164 in subroutine find_var'
        print *, '!!! when analysing the level dimension'
        print *, '!!! of the variable ',name_var(n)
      endif
      go to 160 ! skip name
      172 continue ! from the reading of the level values

      
    ! inquiring the values for the time dimension 
    !--------------------------------------------
      ret3 = nf90_inquire_dimension(ncid,dimids(iddim(3)), &
                                    namedim,length)
      ret4 = nf90_inq_varid(ncid,namedim,timeid)
      if (ret3 == NF90_NOERR .and. ret4 == NF90_NOERR) then
        if (allocated(valuestime)) deallocate(valuestime)
        allocate (valuestime(length),stat=ret5)
        if (ret5 == 0) then
          ret6 = nf90_get_var(ncid,timeid,valuestime)
          if (ret6 == NF90_NOERR) then
            do t = 1,length
              if(valuestime(t) == ts) then
                go to 173 ! good analysis, t is the timestep required
              endif
            enddo
          endif
        endif
      endif
              
      if ( verb .ge. 1 ) then
        print *, ' '
        print *, '!!! Error 165 in subroutine find_var'
        print *, '!!! when analysing the level dimension'
        print *, '!!! of the variable ',name_var(n)
      endif
      go to 160 ! skip name
      173 continue ! from the reading of the level values



    ! retriving the data 
    !--------------------
      if (allocated(values)) deallocate(values)
      if (iddim(1) > iddim(2)) then
        allocate (values(jmax,imax),stat=ret3)
      else
        allocate (values(imax,jmax),stat=ret3)
      endif
      if (ret3 /= 0) then
        if (verb .ge. 1) then
          print *, ""
          print *, '!!! Error 159 in subroutine find_var'
          print *, '!!! when allocating values(jmax,imax), program stop'
        endif
        STOP 159
      endif
      
      starts(iddim(1)) = 1
      starts(iddim(2)) = 1
      starts(iddim(3)) = t
                                    
      counts(iddim(1)) = imax
      counts(iddim(2)) = jmax
      counts(iddim(3)) = 1
                    
      if (ndims == 4) then
        starts(iddim(4)) = l
        counts(iddim(4)) = 1
                      
        ret4 = nf90_get_var(ncid,varid,values,starts,counts)
      else 
        ret4 = nf90_get_var(ncid,varid,values,starts(1:3),counts(1:3))
      endif
                    
      if (ret4 /= NF90_NOERR) then
        if (verb .ge. 1) then
          print *, ""
          print *, '!!! Error 166 in subroutine find_var'
          print *, '!!! The values for variable ', name_var(n)
          print *, "!!! can't be retrived and is skipped"
        endif
        go to 160 ! skip name
      endif
          
          
    ! checking if there is any missing values in the array
    !------------------------------------------------------
      ret5 = nf90_get_att(ncid,varid,'missing_value',missingValue)
      if (ret5 == NF90_NOERR) then
        if (any(values == missingValue)) then
          if (verb .ge. 1) then
            print *, ""
            print *, '!!! Error 167 in subroutine find_var'
            print *, '!!! The variable ', name_var(n)
            print *, '!!! has missing values and is skipped'
          endif
          go to 160 ! skip name
        endif
      else
        if (verb .ge. 1) then
          print *, ""
          print *, 'Warnings : in subroutine find_var'
          print *, 'missing value attribute not found'
          print *, 'we suppose there is no missing data'
          print *, 'for the variable ', name_var(n)
        endif
      endif


    ! changing the order of dimensions if needed and filling data_var
    !----------------------------------------------------------------
      if (iddim(1) > iddim(2)) then
        if (need_to_flip_lats) then
          if (need_to_flip_lons) then
            data_var(:,:,n) = transpose(values(jmax:1:-1,imax:1:-1))
          else
            data_var(:,:,n) = transpose(values(jmax:1:-1,:))
          endif
        else
          if (need_to_flip_lons) then
            data_var(:,:,n) = transpose(values(:,imax:1:-1))
          else
            data_var(:,:,n) = transpose(values(:,:))
          endif
        endif
      else
        if (need_to_flip_lats) then
          if (need_to_flip_lons) then
            data_var(:,:,n) = values(imax:1:-1,jmax:1:-1)
          else
            data_var(:,:,n) = values(:,jmax:1:-1)
          endif
        else
          if (need_to_flip_lons) then
            data_var(:,:,n) = values(imax:1:-1,:)
          else
            data_var(:,:,n) = values(:,:)
          endif
        endif
      endif


      readflag(n) = .TRUE.
      if ( verb .ge. 2 ) then
        print *,''
        print *, 'The variable ', name_var(n)
        print *, 'has been properly loaded'
      endif
      160 continue
      
    enddo nameloop
  
  end subroutine

end module netcdf_datafile




!#################################################################################
!################## Module to read a datafile using grig format ##################
!#################################################################################
!---------------------------------------------------------------------------------

module grib_datafile



  contains
  
  subroutine open_grib_file
  end subroutine
  
  
  subroutine getgridinfo_grib
  end subroutine  
  
  
  subroutine getdata_grib
  end subroutine
    

end module grib_datafile





!#################################################################################
!########### Module distributing the functions which read the datafile ###########
!#################################################################################
!---------------------------------------------------------------------------------

module read_datafile
! Abstract: the only purpose of this module is to simplify the coding
! when modules which read new formats are added

  use netcdf_datafile
  use namelist_trk

  contains
  
  subroutine open_file
  
    implicit none
    
    if (fileinfo%file_type == "netcdf") then
      call open_netcdf_file
    endif
    
  end subroutine
  
  
  subroutine getgridinfo
  
    implicit none
    
    if (fileinfo%file_type == "netcdf") then
      call getgridinfo_netcdf
    endif
    
  end subroutine  
  
  
  subroutine getdata(ts)
  
    implicit none
    integer ts
  
    if (fileinfo%file_type == "netcdf") then
      call getdata_netcdf(ts)
    endif
    
  end subroutine
    

end module read_datafile


!#################################################################################
!######################### Module for the file fort.15 ###########################
!#################################################################################
!---------------------------------------------------------------------------------

module time_step

  use namelist_trk

  implicit none

  integer, save, allocatable :: ts_num(:), ts_tim(:) ! previously ltix and ifhours
  integer, save :: nts !previously ifhmax

  contains

  !---------------------------------------------------------------------------------
  !------------- Subroutine distributive to retrive the timesteps ------------------
  !--------------------------------------------------------------------------------- 
  subroutine readts
  ! ABSTRACT: timesteps could be retrive through different ways:
  ! - a file fort.15 containing the values
  ! - using all the time values contain in the netcdf data file
  
    implicit none
    logical ex
    
    inquire(file = "fort.15", exist = ex)
    
    if(ex) then
      call readfort15
    else if (fileinfo%file_type == "netcdf") then
      call readts_netcdf
    else 
      if (verb .ge. 1) then
        print *, " "
        print *, "!!! ERROR 125 in readts"
        print *, "!!! if the file type is not netcdf, a file fort.15"
        print *, "!!! giving the timestep is needed"
      endif
      STOP 125
    endif
    
    if (verb .ge. 2) then
      print *, "After call readts, the number of time step read is: "
      print *, nts
    endif
    
    
  end subroutine
      
  
  !---------------------------------------------------------------------------------
  !-------------------- Subroutine reading the file fort.15 ------------------------
  !---------------------------------------------------------------------------------    
  subroutine readfort15
  ! previously read_fhours, it reads in a text file (fort.15) that  
  ! contains the forecast times, the format is :
  !    1    0
  !    2    6
  !    3   12
  !  ...  ...
  ! 9999 9999
  ! The first column is the timestep number, the second is the time 
  ! since startdate, given in the units time_unit.
  ! ts_num is filled by the 1st column of fort.15
  ! ts_tim is filled by the 2st column of fort.15
  ! nts is the number of lines (timesteps) in fort.15
  
  
  ! Rq need to create ts_num and ts_tim outside of that

  
    implicit none  
    integer, parameter :: iunit_ts=15 ! fort.15 file
    integer its_num(9999), its_tim(9999)
    integer i, err_a1, err_a2
    
    
    if ( verb .ge. 3 ) then
      print *,'Reading fort.15 file ...'
    endif
          
    nts = 0
    do while (.true.)
      nts = nts + 1
      read (iunit_ts, "(i4,1x,i4)", end=120) its_num(nts), its_tim(nts)      
      if ( verb .ge. 3 ) then
        print "(A,i4,A,i4)", "Read fort.15 : its_num = ", its_num, &
                                            "its_tim = ", its_tim
      endif      
    enddo   
    120 continue


    if (allocated(ts_num))  deallocate (ts_num)
    if (allocated(ts_tim))  deallocate (ts_tim)       
    allocate (ts_num(nts),stat=err_a1)
    allocate (ts_tim(nts),stat=err_a2)
    if (err_a1 /= 0 .or. err_a2 /= 0) then 
      if (verb .ge. 1) then
        print *, ""
        print *, "!!! ERROR 121: allocating ts_num or ts_tim failed"
      endif
      STOP 121
    endif

    
    ts_num = its_num(1:nts)
    ts_tim = its_tim(1:nts)
    
   
    if (minval(ts_tim(2:nts) - ts_tim(1:nts-1)) .le. 0) then
      if ( verb .ge. 1 ) then
        i = minloc(ts_tim(2:nts) - ts_tim(1:nts-1), dim = 1) + 1
        print *,' '
        print *, "!!! ERROR 122 when reading fort.15 :"
        print "(A,i4)", "!!! the time read at line", i
        print *, "!!! is not greater than the previous one :"
        print "(A,i4)", "!!! its_tim(i)=   ", its_tim(i)
        print "(A,i4)", "!!! its_tim(i-1)= ", its_tim(i-1)
        print *,'!!! STOPPING EXECUTION'
      endif
      STOP 122
    endif  
    
  end subroutine
  
  
  
  !---------------------------------------------------------------------------------
  !-------------- Subroutine reading timesteps from a netcdf datafile --------------
  !---------------------------------------------------------------------------------  
  subroutine readts_netcdf
  ! This subroutine is called only for netcdf datafiles and
  ! if a fort.15 file is not provided by the user
  ! It reads all the timesteps in the datafile up to 9999
  ! and complete ts_num, ts_tim and nts
  
  ! we suppose that startdate is the startdate of the datafile
  ! and that the unit is time_units

    USE netcdf
    
    implicit none

    integer ncid, dimid, length, i, time, pos
    integer, allocatable :: values(:)
    integer err_a1, err_a2, ret1, ret2, ret3, ret4, ret5
    character(len = 35) units, startdate, c, buffer
    character(len = 7) time_unit
    integer :: Y,M,D,h,n,s = 0


!   opening and reading the netcdf data file
    ret1 = nf90_open('fort.11', nf90_nowrite, ncid)
    if (ret1 /= nf90_noerr) then
      print *,''
      print *,'!!! ERROR 111 : opening fort.11 failed with netcdf tools'
      STOP 111
    endif

    ret1 = nf90_inq_dimid(ncid, 'time', dimid)
    if (ret1 == nf90_noerr) then
      ret2 = nf90_inquire_dimension(ncid, dimid, len = length)
      if (ret2 == nf90_noerr) then
        allocate(values(length), stat = ret3)
        if (ret3 == 0) then
          ret4 = nf90_get_var(ncid, dimid, values)
          if (ret4 == nf90_noerr) then
            ret5 = nf90_get_att(ncid,dimid,'units', units)
            if (ret5 /= nf90_noerr .and. verb .ge. 2) then
              print *, ''
              print *, "WARNING : units can't be check in the datafile"
            endif
            go to 101
          endif
        endif
      endif
    endif

    print *,''
    print *,'!!! ERROR 112 reading time dimension in fort.11 failed'
    STOP 112

    101  continue

    
    ! Test if the time units in the netcdf file correspond to the requirements
    ! "units" string has the format: "time_units since startdate"
    if (ret5 == nf90_noerr) then
      pos = scan(units, " ")
      time_unit = units(1:pos-1)
      buffer = units(pos+1:)
      
      if (time_unit /= fileinfo%time_unit) then
        print *, ""
        print *, "!!! ERROR 113: The time unit in the netcdf file doesn't"
        print *, "!!! correspond to the one provided in the namelist"
        STOP 113
      endif
      
      pos = scan(buffer, " ")
      c = buffer(1:pos-1)
      buffer = buffer(pos+1:)
      
      !reading the date
      if (c .eq. "since") then

        pos = scan(buffer, "-")
        read(buffer(1:pos-1), "(I4)") Y
        buffer = buffer(pos+1:)
        
        pos = scan(buffer, "-")
        read(buffer(1:pos-1), "(I2)") M
        buffer = buffer(pos+1:) 

        pos = scan(buffer, " ")
        read(buffer(1:pos-1), "(I2)") D
        buffer = buffer(pos+1:) 
        
        pos = scan(buffer, ":")
        read(buffer(1:pos-1), "(I2)") h
        buffer = buffer(pos+1:) 

        pos = scan(buffer, ":")
        read(buffer(1:pos-1), "(I2)") n
        read(buffer(pos+1:), "(I2)") s
        
        write(startdate, "(I0.4,A,I0.2,A,I0.2,A,I0.2,A,I0.2,A,I0.2)") &
                  Y,"-", M,"-", D," ", h,":", n,":", s
                   
        if (startdate /= trkrinfo%startdate) then
          print *, ""
          print *, "!!! ERROR 114: The startdate in the netcdf file"
          print *, "!!! in the namelist"   
          STOP 114
        endif
      
      endif
    endif 


    if (length > 9999 .and. verb .ge. 2) then
      print *, ""
      print *, "WARNING: There are too many timesteps in the datafile"
      print *, "only the 9999 first of ", length, " will be considered"
    endif
    
    nts = min(length, 9999)
    
    if (allocated(ts_num))  deallocate (ts_num)
    if (allocated(ts_tim))  deallocate (ts_tim) 
    allocate (ts_num(nts),stat=err_a1)
    allocate (ts_tim(nts),stat=err_a2)
    if (err_a1 /= 0 .or. err_a2 /= 0) then 
      if (verb .ge. 1) then
        print *, ""
        print *, "!!! ERROR 116: allocating ts_num or ts_tim failed"
      endif 
      STOP 116
    endif
    
    
    do i = 1,nts
      if (abs(values(i)) < 9999) then
        ts_num(i) = i
        ts_tim(i) = values(i)
      else
        print *,''
        print *,'!!! ERROR 117 when filling fort.15'
        print *,'!!! one of the time value exceed the format'
        STOP 117
      endif
    enddo
    
  
  end subroutine



 
end module time_step
  

















