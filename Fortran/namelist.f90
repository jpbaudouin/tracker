!#######################################################################
!### DEFINITION OF THE MODULES FOR THE READING OF THE NAMELIST FILE ####
!#######################################################################


!25-10-2017 Creation
!29-11-2017 Working version
!12-12-2017 1st Commit

! RQ for further developments:
!   Are ERR/END needed in read statement? -> the program stops if an error occured and err/end not specified
!   Create a subroutine writing the namelist ?? --> do it rather with R
!   write test / print_namelist for namelist_tc module
!   The boundaries domain can be -999.0
!   To specify boundaries spanning accross the Greenwitch meridian: trkrinfo%westbd < 0

! CONTAINS :
!   The modules which define extended namelists:
!     - "namelist_tc", for tropical cyclones
!         It includes a namelist "tc_type" and the subroutine read_nlists_tc
! 
!   The module "namelist_trk" for the reading of general parameters :
!     It includes the namelists : "trackparam", "datafile", 
!                                 "fields_names", "verbose"
!                 the general subroutine "read_nlists"
!                 the subroutine print_namelist




!###########################################################################################
!########################## Modules for extended namelists #################################
!###########################################################################################
!-----------------------------------------------------------------------

module namelist_tc

  implicit none

  ! Type declaration ---------------------------------------------------
  !---------------------------------------------------------------------
  type tc_type
    real        mslpthresh ! min mslp gradient to be maintained
    real        v850thresh ! min avg 850 Vt to be maintained
    character*1 haseflag   ! Will phase be determined
                                            ! (y/n)                   
    character*4 hasescheme ! What scheme to use:
                                    ! cps  = Hart's cyclone phase space
                                    ! vtt  = Vitart                    
                                    ! both = Both cps and vtt are used
    real        wcore_depth  ! The contour interval (in deg K)
                             ! used in determining if a closed contour 
                             ! exists in the 300-500 mb T data, for
                             ! use with the vtt scheme.
    character*1 structflag  ! Will structure be
                                            ! analyzed (y/n)?
    character*1 ikeflag     ! Will IKE & SDP be
                                            ! computed (y/n)?
  end type tc_type

  
  ! variables declaration ----------------------------------------------
  !---------------------------------------------------------------------
  type (tc_type)  tc


  contains

  !----- Define now the subroutine which will read the namelist file ---
  !---------------------------------------------------------------------
  subroutine  read_nlists_tc(u_nl, verb)

    implicit none
    
    integer, intent(in)  :: u_nl, verb ! the unit linking to the namelist file
    
    namelist/tc_type/tc
    
    read (u_nl, nml=tc_type)

    if ( verb .ge. 3 ) then
      print *,' '
      print *, 'tc_type :'
      print *, 'mslpthresh :', tc%mslpthresh
      print *, 'v850thresh :', tc%v850thresh
      print *, 'haseflag :', tc%haseflag
      print *, 'hasescheme :', tc%hasescheme
      print *, 'wcore_depth :', tc%wcore_depth
      print *, 'structflag :', tc%structflag
      print *, 'ikeflag :', tc%ikeflag
    endif
                                            
  end subroutine      
end module namelist_tc




!###########################################################################################
!#################################### Main module ##########################################
!###########################################################################################
!-----------------------------------------------------------------------

module namelist_trk

  ! add here the modules defining extended namelist
  use namelist_tc;

  implicit none

  ! Type declaration ---------------------------------------------------
  !---------------------------------------------------------------------
  type trackparam  ! Define a new type for various tracker parmeters
    sequence
    
    character*7   type  ! 'tc', 'midlat', 'WD' or 'other' (more if needed)
    logical       genesis_flag  ! Flag for whether consider the formation of new low
    character*3   detec_minmax  ! tell if looking for a 'min' or a 'max' in the detection field
    real          detec_thresh  ! previously mslpthresh2, min mslp to be detected
    real          contint  ! Contour interval to be used for 
                           ! "midlat" or "tcgen" cases.  ?? not for tracking only   
                           ! what about different fixcenter field ?
                           
    character*19  startdate !previously bcc, byy ... format : YYYY-MM-DD hh:mm:ss
                  ! to be used as origin of fort.15
                      ! (unit define by fort.15, depending on model_start_date)
    integer       output_freq  ! previsouly atcffreq, frequency of output:
                               ! 1 -> every timestep, 2 -> every two timesteps ...
                           
    real          westbd  ! Western boundary of search area
    real          eastbd  ! Eastern boundary of search area
    real          northbd ! Northern boundary of search area
    real          southbd ! Southern boundary of search area

!??  logical     want_oci ! calculate outermost closed isobar
!??  character*1 out_vit  ! Flag for whether to write out vitals, always yes ??

  end type trackparam
  
  !---
  
  type datafile ! Define a new type for various datafile characteristics
    sequence

    character*6   file_type ! previously type, 'netcdf' or 'grib'
    character*8   grid_type ! previously grid, 'global' or 'regional'; merge with modtyp    
    character*7   time_unit ! previously lt_units ('hours' or 'minutes'
                    ! to indicate the units of lead times in grib files)
                    ! now 'seconds', 'minutes', 'hours', 'days', 'months', 'years'.
                    ! and also used to caracterise the unit of fort.15

    character*4   model_name ! previously atcfname (ATCF Name of model (GFSO for GFS...), gmodname ?)
    ! also used for specific parametrisation of hurricane tracking
!??  integer     model_id  ! prevously 'model', integer identifier for model data used (used in get_data_grib only, how useful?) 
    character*19  model_startdate ! previously atcfymdh  used in output : YYYY-MM-DD hh:mm:ss

  end type datafile
  
  !---

  type fields_length
    sequence
    
    integer fixcenter
    integer advection
    integer intensity
    
  end type fields_length

  !---

  type fields_names
    sequence
    
    character*16 detection
    integer     detection_lev
        
    character*16 :: fixcenter(20)
    integer     :: fixcenter_lev(20)
    
    character*16 :: adv_1stdim(20)
    character*16 :: adv_2nddim(20)
    integer     :: adv_lev(20)
    integer     :: adv_factor(20)
    
    character*16 :: intensity(20)
    integer     :: intensity_lev(20)

  end type fields_names
  
  
  ! variables declaration ----------------------------------------------
  !---------------------------------------------------------------------

  type (trackparam)     trkrinfo
  type (datafile)       fileinfo
  type (fields_length)  numfield
  type (fields_names)   fname
  

  ! verbose will be much more limited
  integer, save :: verb            ! Level of detail printed to terminal
                                         ! 0 = No output
                                         ! 1 = Error messages only
                                         ! 2 = Main messages
                                         ! 3 = All 
                                   

  contains

  !---------------------------------------------------------------------------------
  !----------- This subroutine print the values read in the namelist file ----------
  !---------------------------------------------------------------------------------
  subroutine print_namelist

    print *,' '
    print *,'General information read from the namelist :'
    print *, ' '
    print *, 'trackparam'
    print "(A15, A, A)",      ' type', ' : ', &
                                              trkrinfo%type
    print "(A15, A, L1)",     ' genesis_flag', ' : ', &
                                              trkrinfo%genesis_flag
    print "(A15, A, A)",      ' detec_minmax', ' : ', &
                                              trkrinfo%detec_minmax
    print "(A15, A, F20.10)", ' detec_thresh', ' : ', &
                                              trkrinfo%detec_thresh
    print "(A15, A, F20.10)", ' contint',      ' : ', &
                                              trkrinfo%contint
    print "(A15, A, A)",      ' startdate',    ' : ', &
                                              trkrinfo%startdate
    print "(A15, A, I2)",     ' output_freq',  ' : ', &
                                              trkrinfo%output_freq
    print "(A15, A, F8.3)",   ' westbd',       ' : ', &
                                              trkrinfo%westbd
    print "(A15, A, F8.3)",   ' eastbd',       ' : ', &
                                              trkrinfo%eastbd
    print "(A15, A, F8.3)",   ' northbd',      ' : ', &
                                              trkrinfo%northbd
    print "(A15, A, F8.3)",   ' southbd',      ' : ', & 
                                              trkrinfo%southbd
      
    print *, ' '    
        
    print *, 'datafile'
    print "(A15, A, A)", ' file_type',       ' : ', &
                                            fileinfo%file_type
    print "(A15, A, A)", ' grid_type',       ' : ', &
                                            fileinfo%grid_type
    print "(A15, A, A)", ' time_unit',      ' : ', &
                                            fileinfo%time_unit
    print "(A15, A, A)", ' model_name',      ' : ', &
                                            fileinfo%model_name
    print "(A15, A, A)", ' model_startdate', ' : ', &
                                            fileinfo%model_startdate
                                            
    print *, ' '
      
    print *, 'fields_length'
    print "(A15, A, I1)", ' fixcenter', ' :', numfield%fixcenter
    print "(A15, A, I1)", ' advection', ' :', numfield%advection
    print "(A15, A, I1)", ' intensity', ' :', numfield%intensity
      
    print *, ' '
      
    print *, 'fields_names'
    print "(A15, A, A10)", ' detection',       ' : ' , &
                              fname%detection
    print "(A15, A, I10)", ' detection_lev',   ' : ', &
                              fname%detection_lev
    print "(A15, A, 20A10)", ' fixcenter',     ' : ', &
                              fname%fixcenter(1:numfield%fixcenter)
    print "(A15, A, 20I10)", ' fixcenter_lev', ' : ', &
                              fname%fixcenter_lev(1:numfield%fixcenter)
    print "(A15, A, 20A10)", ' adv_1stdim',    ' : ', &
                              fname%adv_1stdim(1:numfield%advection)
    print "(A15, A, 20A10)", ' adv_2nddim',    ' : ', &
                              fname%adv_2nddim(1:numfield%advection)
    print "(A15, A, 20I10)", ' adv_lev',       ' : ', &
                              fname%adv_lev(1:numfield%advection)
    print "(A15, A, 20I10)", ' adv_factor',    ' : ', &
                              fname%adv_factor(1:numfield%advection)
    print "(A15, A, 20A10)", ' intensity',     ' : ', &
                              fname%intensity(1:numfield%intensity)
    print "(A15, A, 20I10)", ' intensity_lev', ' : ', &
                              fname%intensity_lev(1:numfield%intensity)
                              
    print *, ' '
      
    print "(A15, A, I1)", ' verbose ', ' : ', verb
      
  end subroutine


  !---------------------------------------------------------------------------------
  !----------- Define now the subroutine which will read the namelist file ---------
  !---------------------------------------------------------------------------------
  subroutine read_nlists

    ! ABSTRACT: This subroutine simply reads in the namelists that are
    ! provided by the user in the file "namelist"
    implicit none

    integer Y,M,D,h,n,s
    integer iread
    character*5 c
    logical :: err_nl = .FALSE.
        
    namelist/trackparam/trkrinfo
    namelist/datafile/fileinfo
    namelist/fields_length/numfield
    namelist/fields_names/fname
    namelist/verbose/verb

    open(unit = 8, file = "namelist")
    
    read (8,NML=trackparam)
    read (8,NML=datafile) 
    read (8,NML=fields_length)
    read (8,NML=fields_names)
    read (8,NML=verbose) 
    
    
    ! ------ Test if the values are correct :
    !-------------------------------------------------------------------
    print *, ''
    read(trkrinfo%startdate, "(I4,A,I2,A,I2,A,I2,A,I2,A,I2)", &
      iostat = iread) Y, c(1:1), M, c(2:2), D, c(3:3), &
                   h, c(4:4), n, c(5:5), s
    if (iread .eq. 0) then
      if(c /= "-- ::")  then
      err_nl = .TRUE.
      print * ,'trkrinfo%startdate'
      endif
    else
      err_nl = .TRUE.
      print * ,'trkrinfo%startdate'
    endif

    read(fileinfo%model_startdate, "(I4,A,I2,A,I2,A,I2,A,I2,A,I2)", &
      iostat = iread) Y, c(1:1), M, c(2:2), D, c(3:3), &
                   h, c(4:4), n, c(5:5), s
    if (iread .eq. 0) then
      if(c /= "-- ::")  then
      err_nl = .TRUE.
      print * ,'fileinfo%model_startdate'
      endif
    else
      err_nl = .TRUE.
      print * ,'fileinfo%model_startdate'
    endif
    
    if ((fileinfo%file_type /= 'netcdf') .and. &
        (fileinfo%file_type /= 'grib')) then
      err_nl = .TRUE.
      print *, ''
    endif
    
    if ((fileinfo%grid_type /= 'global') .and. &
        (fileinfo%grid_type /= 'regional')) then
      err_nl = .TRUE.
      print *, ''
    endif

    if ((fileinfo%time_unit /= 'seconds') .and. &
        (fileinfo%time_unit /= 'minutes') .and. &
        (fileinfo%time_unit /= 'hours') .and. &
        (fileinfo%time_unit /= 'days') .and. &
        (fileinfo%time_unit /= 'months') .and. &
        (fileinfo%time_unit /= 'years')) then
      err_nl = .TRUE.
      print *, ''
    endif

 
    if ((trkrinfo%detec_minmax /= 'min') .and. &
        (trkrinfo%detec_minmax /= 'max')) then
      err_nl = .TRUE.
      print *, 'trkrinfo%detec_minmax'
    endif
    
    if (trkrinfo%output_freq .le. 0) then
      err_nl = .TRUE.
      print *, 'trkrinfo%output_freq'
    endif

    if (trkrinfo%northbd .gt. 90.) then
      err_nl = .TRUE.
      print *, 'trkrinfo%northbd'
    endif

    if ((trkrinfo%southbd .lt. -90.) .and. &
        (trkrinfo%southbd .ge. -998.)) then
      err_nl = .TRUE.
      print *, 'trkrinfo%southbd'
    endif

    if (trkrinfo%southbd .gt. trkrinfo%northbd) then
      err_nl = .TRUE.
      print *, 'trkrinfo%southbd > trkrinfo%northbd'
    endif
    
    if ((trkrinfo%westbd .lt. -360.) .and. &
       (trkrinfo%westbd .ge. -998.)) then
      err_nl = .TRUE.
      print *, 'trkrinfo%westbd'
    endif

    if (trkrinfo%eastbd .gt. 360. .or. (trkrinfo%eastbd .lt. 0 .and. &
       trkrinfo%eastbd .ge. -998.)) then
      err_nl = .TRUE.
      print *, 'trkrinfo%eastbd'
    endif
    
    if (trkrinfo%eastbd - trkrinfo%westbd .gt. 360.) then
      err_nl = .TRUE.
      print *, 'trkrinfo%eastbd - trkrinfo%westbd > 360'
    endif    

    if (trkrinfo%westbd .gt. trkrinfo%eastbd) then
      err_nl = .TRUE.
      print *, 'trkrinfo%westbd > trkrinfo%eastbd'
    endif
        
    if ((numfield%fixcenter .lt. 0) .or. &
        (numfield%fixcenter .gt. 20)) then
       err_nl = .TRUE.
      print *, 'numfield%fixcenter'
    endif 

    if ((numfield%advection .le. 0) .or. &
        (numfield%advection .gt. 20)) then
       err_nl = .TRUE.
      print *, 'numfield%advection'
    endif

    if ((numfield%intensity .lt. 0) .or. &
        (numfield%intensity .gt. 20)) then
       err_nl = .TRUE.
      print *, 'numfield%intensity'
    endif

    
    if (err_nl) then
      print *, '!!! ERROR 101 in reading namelist :'
      print *, 'The variables above are wrong, see what has been read:'  
      call print_namelist
      STOP 101    
    endif

    if (verb .ge. 3) then
      call print_namelist
    endif

  
    !-------------------------------------------------------------------
    ! for extended namelist, add the call here
    ! it should depend on trkrinfo%type
    !-------------------------------------------------------------------    
    if (trkrinfo%type == 'tc') then
      call read_nlists_tc(8, verb)
    endif
    
    !-------------------------------------------------------------------
    
    close(8)
    
  end subroutine
  
end module namelist_trk

