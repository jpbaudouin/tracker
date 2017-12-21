!#######################################################################
!########## DEFINITION OF THE CORE PARAMETERS AND SUBROUTINES ##########
!#######################################################################


!25-10-2017 Creation
!12-12-2017 1st Commit -> pi, dtr
!14-12-2017 2nd commit -> check_closed_contour, first_ges_center
!21-12-2017 3rd commit -> module no_dependance, get_neighbours
!                         get_ij_bounds, find_minmax, fixcenter
!                         correction check_contour (check_closed_contour), first_ges_center

! Is there a better way to change max_object than recompile it?
?? ! Watch out: mask_detec and masked_out have opposite meaning for their values
?? ! all .true. are considered for mask_detec, while all .false. are for masked_out


! CONTAINS :
!   The module no_dependance which includes subroutines without dependance 
!   to parameters from other modules:
!   - "get_neighbours", which returns a list of neighboring points around
!       a central one

!   The module core_tracker with different parameters and the subroutines:
!   - "get_ij_bounds", which reduces the array size for local analysis,
!       by returning the indeces of a subgrid which contains all the 
!       points within a certain distance of a central one.




!###########################################################################################
!##################### Module including subrtoutine without dependances ####################
!###########################################################################################
!-----------------------------------------------------------------------------------

module no_dependance

  contains
  

  !---------------------------------------------------------------------------------
  !-------------- Subroutine returning a list of neighboring points ----------------
  !---------------------------------------------------------------------------------
  subroutine get_neighbours(pt_i, pt_j, grid, imax, jmax ptn, ptn_i, ptn_j)
  
  ! This function return a list of ptn points represented by their 
  !   indices (ptn_i, ptn_j) that are around a point of index (pt_i, pt_j).
  ! It considers at maximum the 8 neighboring points of (pt_i, pt_j).
  ! If the grid is "global", it considers a possble wrapping along i
  ! If the grid is "regional", and (pt_i, pt_j) is on its edge,
  !   it will reduce ptn to the number of existing neighbours (5 or 3)
  ! There is no dependance on the original grid, in order to be used on 
  !   local grid as given by get_ij_bounds

  
  ! INPUT
  ! pt_i      Index along i of the central point (0 < pt_i <= i_max)
  ! pt_j      Index along j of the central point (0 < pt_j <= j_max)
  ! grid      the type of grid for possible wrapping, either "regional"
  !             or "global"
  ! imax      The maximum number of point along i on the grid
  ! jmax      The maximum number of point along j on the grid

  ! OUTPUT
  ! ptn       number of neighboring points (8, 5 or 3)
  ! ptn_i     List of index along i of the neighboring points
  ! ptn_j     List of index along j of the neighboring points
  
  
    implicit none
  
    integer, intent(in)  :: pt_i, pt_j
    integer, intent(out) :: ptn, ptn_i(8), ptn_j(8)
  
    integer pt_ip1, pt_im1, pt_jp1, pt_jm1


  !-------------------------------------------------  
  !-- Determine i / j of the 8 neighboring points --
  !-------------------------------------------------
    pt_ip1 = pt_i(ipt) + 1
    pt_im1 = pt_i(ipt) - 1
    pt_jp1 = pt_j(ipt) + 1
    pt_jm1 = pt_j(ipt) - 1
    if (grid == 'global') then ! wrapping
      if (pt_i(ipt) == 1)    pt_im1 = imax
      if (pt_i(ipt) == imax) pt_ip1 = 1
    endif

    ! list of the 8 neighboring points
    ptn_i = [pt_ip1, pt_ip1   , pt_ip1, pt_i(ipt), &
             pt_im1, pt_im1   , pt_im1, pt_i(ipt)  ]
    ptn_j = [pt_jp1, pt_j(ipt), pt_jm1, pt_jm1,    &
             pt_jm1, pt_j(ipt), pt_jp1, pt_jp1     ]
    ptn = 8

  !---------------------------------------------------------------------
  !-- if the border of a grid is hit, only consider the reamining points
  !-- (if model is global, only Sth and Nth or considered) -------------
  !---------------------------------------------------------------------
    ! along j
    if (pt_j(ipt) == 1) then
      ptn = 5
      ptn_i(1:ptn) = [pt_ip1, pt_ip1   , pt_im1   , pt_im1, pt_i(ipt)]
      ptn_j(1:ptn) = [pt_jp1, pt_j(ipt), pt_j(ipt), pt_jp1, pt_jp1   ]
    else if (pt_j(ipt) == jmax) then
      ptn = 5
      ptn_i(1:ptn) = [pt_ip1   , pt_ip1, pt_i(ipt), pt_im1, pt_im1   ]
      ptn_j(1:ptn) = [pt_j(ipt), pt_jm1, pt_jm1   , pt_jm1, pt_j(ipt)]
    endif

    ! along i and corners
    if (grid == 'regional') then
      if (pt_i(ipt) == 1) then
        if (pt_j(ipt) == 1) then
          ptn = 3
          ptn_i(1:ptn) = [1,2,2]
          ptn_j(1:ptn) = [2,2,1]
        else if (pt_j(ipt) == jmax) then
          ptn = 3
          ptn_i(1:ptn) = [1,2,2]
          ptn_j(1:ptn) = [jmax-1,jmax-1,jmax]
        else
          ptn = 5
          ptn_i(1:ptn) = [pt_ip1,pt_ip1   ,pt_ip1,pt_i(ipt),pt_i(ipt)]
          ptn_j(1:ptn) = [pt_jp1,pt_j(ipt),pt_jm1, pt_jm1  ,pt_jp1]
        endif
      else if (pt_i(ipt) == imax) then
        if (pt_j(ipt) == 1) then
          ptn = 3
          ptn_i(1:ptn) = [imax,imax-1,imax-1]
          ptn_j(1:ptn) = [2,2,1]
        else if (pt_j(ipt) == jmax) then
          ptn = 3
          ptn_i(1:ptn) = [imax,imax-1,imax-1]
          ptn_j(1:ptn) = [jmax-1,jmax-1,jmax]
        else
          ptn = 5
          ptn_i(1:ptn) = [pt_i(ipt),pt_im1,pt_im1   ,pt_im1,pt_i(ipt)]
          ptn_j(1:ptn) = [pt_jm1   ,pt_jm1,pt_j(ipt),pt_jp1,pt_jp1]
        endif
      endif
    endif
             
  end subroutine

end module




!###########################################################################################
!#################### Module including core subrtoutines and parameters ####################
!###########################################################################################
!----------------------------------------------------------------------------------

module core_tracker
  ! core subroutine and parameters
  
  use data_param
  use namelist_trk

  implicit none
    
  real, parameter :: pi = 4. * atan(1.)
  real, parameter :: dtr = 4. * atan(1.)/180.0 !pi/180

  integer, parameter :: max_object = 5000
  
?!  masked_out and ilast_object
  
  
  contains
  


  !---------------------------------------------------------------------------------
  !-------------- Subroutine reducing array size for local analysis ----------------
  !---------------------------------------------------------------------------------

  subroutine get_ij_bounds (geslon, geslat, rad, &
               ibeg, jbeg, iend, jend, ret)
  
  ! It returns the indices of a subgrid that includes all the points
  !   within a distance rad of the point of coordinate (geslon, geslat)
  ! Note that the point (geslon, geslat) doesn't need to be on the main grid
  ! Note also that the points at the edge of this subgrid are all  
  !   farther than rad of the central point.
  ! Effectively, ibeg, jbeg, iend, jend represent the index of the first
  !   point farther than rad towards, respectively the West, the North, 
  !   the East and the South of the point (geslon, geslat)
  
  ! INPUT
  ! geslon    longitude coordinate of the central point
  ! geslat    latitude  coordinate of the central point
  ! rad       the radius around the central point
  
  ! OUTPUT
  ! ibeg      first index along i of the subgrid (wrapping => ibeg > iend)
  ! jbeg      first index along j of the subgrid
  ! iend      last index along i of the subgrid (wrapping => ibeg > iend)
  ! jend      last index along j of the subgrid
  ! ret       return success of the subroutine
  !           0   : no error
  !           71  : rad is too small and the first point tested is 
  !                   already at a distance superior to rad
  !           72  : rad is likely too large, as all indexes along i are considered
  !           < 0 : In case of a regional grid, the distance between 
  !                   the edge and the central point is inferior to rad 
  !                 In that case, the edge of the subgrid is fixed to
  !                   the edge of the main one

  
    real,    intent(in)  :: geslon, geslat, rad
    integer, intent(out) :: ibeg, jbeg, iend, jend, ret
    integer fix_i, fix_j, ip, jp, i, j
    real dist, degrees

    ret = 0
    
  !-----------------
  !-- to the East --
  !-----------------
    fix_i = int((geslon - glonmin + 2*dx)/dx))
    do ip=fix_i, fix_i + imax - 1
      
      if (ip <= imax) then
        i = ip
      else if (fileinfo%grid_type == 'regional') then
        ret = -1
        iend = imax
        exit
      else ! global model
        if (ip >= fix_i + imax/2 + 1) then
        ! it's about to loop around the earth and to get closer to gespoint
          iend = imax ! ibeg will be 1 in its loop
          ret = 72
          exit
        else ! wrapping
          i = ip - imax
        endif
      endif
      
      call calcdist(geslon,geslat,glon(i),geslat,dist,degrees)
      
      if (dist > rad) then      
        if (i == fix_i) then
          ret = 71
        endif       
        iend = i
        exit
      endif
    enddo

  !-----------------
  !-- to the West --
  !-----------------
    fix_i = int((geslon - glonmin + dx)/dx))
    do ip=fix_i, -(fix_i + imax - 1), -1
      
      if (ip >= 1) then
        i = ip
      else if (fileinfo%grid_type == 'regional') then
        ret = ret - 10
        ibeg = 1
        exit
      else ! global model
        if (ip <= -(fix_i + imax/2 + 1)) then
        ! it's about to loop around the earth and to get closer to gespoint
          ibeg = 1 ! iend will be imax in its loop
          ret = 72
          exit
        else ! wrapping
          i = ip + imax
        endif
      endif
      
      call calcdist(geslon,geslat,glon(i),geslat,dist,degrees)
      
      if (dist > rad) then      
        if (i == fix_i) then
          ret = 71
        endif       
        ibeg = i
        exit
      endif
    enddo        

  !-----------------
  !-- to the South--
  !-----------------
    fix_j = int((geslat - glatmin + 2*dy)/dy))
    do jp=fix_j, fix_j + jmax - 1
      
      if (jp <= jmax) then
        j = jp
      else
        if (fileinfo%grid_type == 'regional') then
          ret = ret - 100
        endif ! this is not a problem in case of a global grid
        jend = jmax
        exit
     endif
      
      call calcdist(geslon,geslat,geslon,glat(j),dist,degrees)
      
      if (dist > rad) then      
        if (j == fix_j) then
          ret = 71
        endif       
        jend = j
        exit
      endif
    enddo  
    
  !-----------------    
  !-- to the North--
  !-----------------
    fix_j = int((geslat - glatmin + dy)/dy))
    do jp=fix_j, -(fix_j + jmax - 1), -1
      
      if (jp >= 1) then
        j = jp
      else
        if (fileinfo%grid_type == 'regional') then
          ret = ret - 1000
        endif ! this is not a problem in case of a global grid
        jbeg = 1
        exit
      endif 
         
      call calcdist(geslon,geslat,geslon,glat(j),dist,degrees)
      
      if (dist > rad) then      
        if (j == fix_j) then
          ret = 71
        endif       
        jbeg = j
        exit
      endif
    enddo  

 
  end subroutine
  
  

  !---------------------------------------------------------------------------------
  !----------- Subroutine searching for a closed contour around a point ------------
  !---------------------------------------------------------------------------------

  subroutine check_contour (ptx_i, ptx_j, fxy, cmaxmin, contour_min, &
                forced, mask_object, mask_object_cc, &
                contour_max, rcontour_max,ret)

! all object should have a minimum deepness of trkrinfo%contint

    ! possible points at the border: either in another object's area or with
    ! a lower value than a point of the object which have a value > contour)
    
    
! masked_out might be an argument if different than "area influenced by other object"
! It seems it's values are updated in the subroutine'
! merge closed_contour with ret (icccret)
! num_requested_contours modified into logical multi_contour (not implemented yet)
! get_last_isobar_flag is supressed


    implicit none

    integer,     intent(in) :: ptx_i, ptx_j
    real,        intent(in) :: fxy(:,:), contour_min
    character*3, intent(in) :: cmaxmin
    logical,     intent(in) :: forced

    logical, intent(out) :: mask_object(:,:)
    logical, intent(out) :: mask_object_cc(:,:)
    real,    intent(out) :: contour_max, rcontour_max
    integer, intent(out) :: ret = 0

    real :: xcentval = fxy(ix,jx)
    real contour, dist, degree
    integer ipt, pt_max, ptn, iptn
    integer pt_i(imax*jmax), pt_j(imax*jmax), ptn_i(8), ptn_j(8)
    logical contour_close, on_border
    integer mm ! 1: max, -1: min


    if ( verb .ge. 3 ) then
      print *,' '
      print *,'*-----------------------------------------------*'
      print *,'* Top of check_closed_contour, ix= ',ix,' jx= ',jx
      print *,'*-----------------------------------------------*'
      print *,'cmaxmin: ', cmaxmin
      print *,'fxy(ix,jx)= ', fxy(ix,jx)
      print *,'minmial contour= ', contour_min
    endif


    ! to simplify the if statements
    if (cmaxmin == 'min') then
      mm = -1
    else
      mm = 1
    endif

  !--------------------
  !-- initialisation --
  !--------------------
    ipt = 1 ! current index
    pt_max = 1 ! maximum index
    pt_i(1) = ptx_i ! the points to be treated
    pt_j(1) = ptx_j 

    mask_object(ptx_i, ptx_j) = .false.

    !nptb = 0 ! number of possible points at the border

    contour_close = .true.


  !-------------------------------------
  !-- Loop on the points (pt_i, pt_j) --
  !-------------------------------------
    scan_contour: do while (ipt <= pt_max)
      if ( verb .ge. 3 ) then
        print *, " "
        print *, " the following point is beeing treated: "
        print *, " -lon: ", glon(pt_i(ipt)), " -lat: ", glat(pt_j(ipt))
      endif 
    
    
      call get_neighbours(pt_i(ipt), pt_j(ipt), fileinfo%grid_type, &
              imax, jmax, ptn, ptn_i, ptn_j)
          
      ! we hit the border of a regional grid
      if (fileinfo%grid_type == 'regional' .and. ptn /= 8) then
      
        if ( verb .ge. 3 ) then
          print *,'We reach the border of a regional grid'
        endif
        
        if ( (contour_min-fxy(pt_i(ipt), pt_j(ipt)))*mm > 0)  then
        
          if ( verb .ge. 3 ) then
            print *,'The minimum contour is not closed'
          endif
          
          contour_close = .false.
          
          if (.not. forced) then
            ret = 61
            exist scan_contour
          endif
          
        else ! the object is deep enough, the point is on the border 
          !nptb = nptb + 1
          !ptb_i(nptb) = pt_i(ipt)
          !ptb_j(nptb) = pt_j(ipt)
          ptn = 0 ! its neighbours are not considered
        endif
        
      endif          


    !-----------------------------------
    !-- scan the 8 neighboring points --
    !-----------------------------------
      scan_neighbour : do iptn=1,ptn
      
      
        if (masked_out(ptn_i(iptn), ptn_j(iptn))) then
          if ( verb .ge. 3 ) then
            print *,'We reach the border of another object'
          endif
          
          if ( (contour_min-fxy(pt_i(ipt), pt_j(ipt)))*mm > 0 ) then
          
            if ( verb .ge. 3 ) then
              print *,'The minimum contour is not closed'
            endif
            
            contour_close = .false.

          if (.not. forced) then
            ret = 62
            exist scan_contour
          endif
            
          else ! the object is deep enough, the point is on the border 
            !nptb = nptb + 1
            !ptb_i(nptb) = ptn_i(iptn)
            !ptb_j(nptb) = ptn_j(iptn)
          endif


        else if ((xcentval-fxy(ptn_i(iptn), ptn_j(iptn)))*mm < 0) then
          if ( verb .ge. 3 ) then
            print *, "We reach a point with a lower/higher value ", &
                      "than the limit"
            print *, "fxy(i,j) :" fxy(ptn_i(iptn), ptn_j(iptn))
            print *, "limit: ", xcentval
            print *, "cminmax: ", cminmax
          endif

          if ( (contour_min-fxy(pt_i(ipt), pt_j(ipt)))*mm > 0 ) then
          
            if ( verb .ge. 3 ) then
              print *,'The minimum contour is not closed'
            endif
            
            contour_close = .false.

          if (.not. forced) then
            ret = 63
            exist scan_contour
          endif
            
          else ! the object is deep enough, the point is on the border 
            !nptb = nptb + 1
            !ptb_i(nptb) = ptn_i(iptn)
            !ptb_j(nptb) = ptn_j(iptn)
          endif
          
        ! There is nothing to do if the point is already in mask_object
        else if (.not. mask_object(ptn_i(iptn), ptn_j(iptn))) then
          if (( (contour_min-fxy(pt_i(ipt), pt_j(ipt)))*mm <= 0 ) .or.
              ( (fxy(pt_i(ipt), pt_j(ipt)) - 
                 fxy(ptn_i(iptn), ptn_j(iptn)))*mm > 0 ) then
        
            mask_object(ptn_i(iptn), ptn_j(iptn)) = .true.
            pt_max = pt_max + 1
            pt_i(pt_max) = ptn_i(iptn)
            pt_j(pt_max) = ptn_j(iptn)
            
          else ! the point is on the border
            !nptb = nptb + 1
            !ptb_i(nptb) = ptn_i(iptn)
            !ptb_j(nptb) = ptn_j(iptn)
          endif
          
        endif
      enddo scan_neighbour

      ipt = ipt + 1
      
    enddo scan_contour



    if (.not. contour_close) then
      mask_object_cc = .false.
      mask_object = .false.
      contour_max = xcentval
      rcontour_max = -999.
      return
    endif
    
  
  !------------------------------------------
  !-- determine the highest closed contour --
  !------------------------------------------
    ! It's the lowest value of the points on the border (cminmax=min)
    contour_max = maxval(fxy, mask = mask_object)
    do i = 1,imax
      do j = 1,jmax
        if(mask_object(i,j)) then
        

          call get_neighbours(i, j, fileinfo%grid_type, &
              imax, jmax, ptn, ptn_i, ptn_j)
          
          on_border = .false.
          do iptn=1,ptn
            if(.not. mask_object(ptn_i(iptn),ptn_j(iptn))) then
              on_border = .true.
            endif
          enddo
          
          if(on_border .and. (contour_max - fxy(i,j))*mm < 0) then
            contour_max = fxy(i,j)
          endif
          
        endif
      enddo
    enddo

    
    ! fill the mask
    mask_object_cc = mask_object
    do i = 1,imax
      do j = 1,jmax
        if(mask_object(i,j) .and. fxy(i,j) >= contour_max) then
          mask_object_cc(i,j) = .false.
        endif
      enddo
    enddo


    ! compute the radius    
    nb = 0
    rcontour_max = 0
    do i = 1,imax
      do j = 1,jmax
        if(mask_object_cc(i,j)) then
        

          call get_neighbours(i, j, fileinfo%grid_type, &
              imax, jmax, ptn, ptn_i, ptn_j)
          
          on_border = .false.
          do iptn=1,ptn
            if(.not. mask_object_cc(ptn_i(iptn),ptn_j(iptn))) then
              on_border = .true.
            endif
          enddo
          
          if(on_border .and. (contour_max - fxy(i,j))*mm < 0) then
            call calcdist(glon(ptx_i), glat(ptx_j), glon(i), glat(j), &
                    dist,degrees)
            nb = nb + 1
            rcontour_max = rcontour_max + dist
          endif
          
        endif
      enddo
    enddo

    rcontour_max = rcontour_max/nb


  end subroutine
  
    

  
  !---------------------------------------------------------------------------------
  !--------------- Subroutine searching for a minmax around a point---- ------------
  !---------------------------------------------------------------------------------
  subroutine find_minmax(geslon, geslat, fxy, minmax, ptc_i, ptc_j, &
                clon, clat, cval)
  
  ! This subroutine looks for either local minimum or maximum of a field 
  !   fxy around a point determine by its coordinate (geslon, geslat).
  ! The first step consists in reducing the area studied to a square
  !   including all points within a certain distance of the point
  !   (cf. get_ij_bounds).
  ! Then, all local minmax within that distance and not influenced by
  !   another object are detected
  ! Finally, the minmax selected is the one maximising the formula 
  !   (val-mean(A)/std(A))*((dist-rad)/rad)**2 where val is its value,
  !   dist its distance to the central point, rad the maximum of this
  !   distance and A the area of search.
  ! If no values are found, clon and clat = -9999.
  
  
?? I am still hesitating if I send back the lat/lon or the indices => it depends if I want to have a better aproximation (the bicubic interpolation) or if I need the centers to be on the grids (ex in fixcenter to include them in the area influenced by the storm). As for now I choose to have both. (but the bicubic interpolation is not yet implemented)

?? ! I can add a bicubic interpolation over the 8 surronding points of the selected minmax to precise its position and its value

! previously find_maxmin     
! there is some special case when cparm = zeta, or cparm = vmag
    ! but in subroutine tracked only zeta, hgt and slp are used
    ! tmp is used in get_vtt_phase
    ! vmag in get_uv_center !!!! there glon, glat, glon/lat-min/max change. Not implemented yet.
    ! and again zeta in get_zeta_values
! pt_i,pt_j is guesslon, guesslat, where to look at
! glat, glon are rlonv, rlatv
! compflag merged with ret (ifmret)
! ptc_i,ptc_j is ctlon, ctlat, the center found, and xval the value
! cmodel_type is trkrinfo%gridtype
! ret is ifmret

?? !rad and minmax in namelist

    real, intent(in) :: geslon, geslat, fxy(:,:)
    character*3, intent(in) :: minmax
    
!    real, intent(out) :: clon, clat
    integet, intent(out) :: ptc_i, ptc_j
    real, intent(out) :: cval
    
    integer ibeg, jbeg, iend, jend ! the boundaries of the local grid
    real, allocatable :: fxy2 ! local version of fxy
    logical , allocatable :: mask_s ! area of search
    integer iend2, imax2, jmax2, iw, ip ! for wrapping issue
    real dist,degrees, mean, sd
    integer nptm, ptm(2), iptm ! for the possible minmax
    integer, allocatable ptm_i(:), ptm_j(:) ! for the possible minmax
    integer ptm, ptc_i, ptc_j ! effective minmax
    integer ptn, ptn_i(8), ptn_j(8), iptn ! the neighbours
    logical local_minmax
    
    
    integer ret
  
    ptm ! possible minmax
    ptm_i, ptm_j ! effective minmax



  !------------------------------------------------
  !-- Set up the mask where to look for a minmax --
  !-- It reduces the area studied -----------------
  !------------------------------------------------
    call get_ij_bounds(geslon, geslat, rad, ibeg, jbeg, iend, jend, ret)

    if (ret == 71) then
      if ( verb .ge. 1 ) then
        print *, ""
        print *, "!!! ERROR 71 in find_minmax, from get_ij_bounds:"
        print *, "!!! The caracteristic size of the object provided by"
        print *, "!!! the user, rad (", rad, ")"
        print *, "!!! is smaller than a gird point."
        print *, "!!! dx: ", dx, " dy: ", dy
        print *, "!!! No minmax can be found, tracking ends"
      endif
      STOP 71
    endif

    if (ret == 72) then
      if ( verb .ge. 1 ) then
        print *, ""
        print *, "!!! ERROR 72 in find_minmax, from get_ij_bounds:"
        print *, "!!! The caracteristic size of the object provided by"
        print *, "!!! the user, rad (", rad, ")"
        print *, "!!! is larger than 1/2 circonference of the Earth"
        print *, "!!! Too many minmax can be found, tracking ends"
      endif
      STOP 72
    endif
          
    if ( verb .ge. 2 .and. ret < 0) then
      print *, ""
      print *, " WARNING in find_minmax, from get_ij_bounds:"
      print *, " The edge of a regional grid has been hit"
      print *, " This might lead to some problem in finding a minmax"
    endif
    
        
    if (ibeg > iend) then ! wrapping, Rq, ibeg > 0 by definition
    
      allocate(mask_s(iend-ibeg+1+imax, jend-jbeg+1))
      allocate(fxy2(iend-ibeg+1+imax, jend-jbeg+1))
      
      fxy2((imax-ibeg+2):(iend-ibeg+1+imax),:) = fxy(1:iend, jbeg:jend)
      fxy2(1:(imax-ibeg+1),:) = fxy(ibeg:imax, jbeg:jend)
      
      iend2 = iend + imax ! for the next loop
      imax2 = iend-ibeg+1+imax ! for the call to get_neighbours
      jmax2 = jend-jbeg+1 ! for the call to get_neighbours
      
    else
    
      allocate(mask_s(iend-ibeg+1, jend-jbeg+1))
      allocate(fxy2(iend-ibeg+1, jend-jbeg+1))
      
      fxy2 = fxy(ibeg:iend, jbeg:jend)

      iend2 = iend ! for the next loop
      imax2 = iend-ibeg+1 ! for the call to get_neighbours
      jmax2 = jend-jbeg+1 ! for the call to get_neighbours
      
    endif

    
    ! fill mask_s, the mask of search, only trues will be analysed
    mask_s = .not. masked_out
    do i=ibeg,iend2 ! Rq: 0 < ibeg < iend2, wrapping or not
    
      !wrapping
      iw = i
      if (i > imax) iw = i - imax
      
      do j=jbeg,jend
      
        call calcdist(geslon,geslat,glon(iw),glat(j),dist,degrees)
        
        if (dist > rad) then
          mask_s(i - ibeg + 1, j - jbeg + 1) = .false.
        endif
      enddo
    enddo
    
    ! To compute the weights at the last step before mask_s changes
    mean = sum(fxy2, mask = mask_s)/count(mask_s)
    sd   = sqrt(sum((fxy2-mean)**2, mask = mask_s)/(count(mask_s)-1))
    

  !-----------------------------
  !-- search for local minmax --    
  !-----------------------------
    nptm = 0 ! number of minmax found
    allocate(ptm_i(imax2*jmax2))
    allocate(ptm_j(imax2*jmax2))
    
    do while (any(mask_s))
    
      if (minmax == 'min') then
        ptm = minloc(fxy2, mask = mask_s)
      else 
        ptm = maxloc(fxy2, mask = mask_s)
      endif
      
      ! It's possible to hit the edge of the grid for some weird case
      call get_neighbours(ptm(1), ptm(2),'regional', imax2, jmax2, &
              ptn, ptn_i, ptn_j)
      
      ! If the edge of the local grid is hit, it's because it has hit
      !   the edge of the main one
      ! If the main grid is regional, it can't be a minmax
      ! If the main grid is global then it should be the northern or
      !   southern border and it's a possible minmax
      if (ptn < 8 .and. fileinfo%grid_type = 'regional') then
        local_minmax = .false.
      else
        local_minmax = .true.      
        do iptn=1,ptn !scan neighbours
        
          if (masked_out(ptn_i(iptn), ptn_j(iptn))) then
            local_minmax = .false.
            exit
          endif
          
          if (minmax == 'min') then
            if (fxy2(ptn_i(iptn), ptn_j(iptn)) < &
                fxy(ptm(1), ptm(2))) then
              local_minmax = .false.
              exit
            endif
          else 
            if (fxy2(ptn_i(iptn), ptn_j(iptn)) > &
                fxy(ptm(1), ptm(2))) then
              local_minmax = .false.
              exit
            endif
          endif
          
        enddo
      endif
      
      mask_s(ptm(1), ptm(2)) = .false.
      
      if (local_minmax) then
        nptm = nptm + 1
        ptm_i(nptm) = ptm(1)
        ptm_j(nptm) = ptm(2)
        
        do iptn=1,ptn
          mask_s(ptn_i(iptn), ptn_j(iptn)) = .false.
        enddo
      endif
      
    enddo  


  !---------------------------------------------------
  !-- Determine the most likely minmax local minmax --    
  !---------------------------------------------------
    if (nptm > 0) then
      do iptm=1,nptm
        call calcdist(geslon,geslat,glon(ptn_i(iptm)),glat(ptn_j(iptm)), &
               dist,degrees)
      
        val(iptm) = ((fxy2(ptm_i(iptm), ptm_j(iptm))-mean)/std)* &
                    ((dist-rad)/rad)**2
      enddo
    
      ptm = maxloc(val)
    
      ptc_i = ptn_i(ptm) + ibeg - 1
      ptc_j = ptn_j(ptm) + jbeg - 1
      
      clon = glon(ptc_i)
      clat = glat(ptc_j)
      cval = fxy2(ptn_i(iptm), ptn_j(iptm))

      if (ptc_i > imax) then  ! wrapping
        ptc_i = ptc_i - imax
      endif
      
    else
      ptc_i = -999
      ptc_j = -999
    
      clon = -9999.
      clat = -9999.
      cval = -9999.
    endif

  
  end subroutine
  


  
  !---------------------------------------------------------------------------------
  !----------------- Subroutine precising the center of an object ------------------
  !---------------------------------------------------------------------------------

  subroutine fixcenter (pt_i, pt_j, mask_object_max, &
                mask_object_min, mask_object, fixlon, fixlat, fixvals, &
                contour_max, rcontour_max)


! add the name of the object in printing


    real, intent(in) :: pt_i, pt_j ! the minmax from field detection
    logical, intent(in) :: mask_object_max(:,:), mask_object_min(:,:)
    
    real, intent(out) :: fixlon,fixlat,fixvals(:)
    logical, intent(out) :: mask_object(:,:)
    
    real, intent(inout) :: contour_max, rcontour_max
    
    real geslon, geslon ! the lon/lat of the center from the field detection
    
      ! the centers in the different fields
    real clon(numfield%fixcenter + 1), clat(numfield%fixcenter + 1)
    integer ptc_i(numfield%fixcenter + 1), ptc_j(numfield%fixcenter + 1)
      
    real fix_lon1, fix_lat1 ! the center after simple average
    real fix_lon, fix_lat ! the center after weighted average
    integer ifc, nfc ! number of fixcenter fields with a minmax
    integer nb, i, j
    real wt, wts, dist, degrees, val
    logical in_mask, on_border, mask_object_notused(imax, jmax)



    if ( verb >= 3 ) then
      print *, " "
      print *, "-------------------------------"
      print *, " At the beginning of fixcenter "
      print *, "-------------------------------"
    endif

     
  !-------------------------------------------------------------
  !-- First, we compute the minmax of the fields in fixcenter --
  !-------------------------------------------------------------
   
    geslon = glon(pt_i)
    geslat = glat(pt_j)
      
    do ifc = 1, numfield%fixcenter
      call findmaxmin(geslon, geslat, data_fixcenter(:,:,ifc), &
              fname%fc_minmax(ifc), ptc_i(ifc), ptc_j(ifc), &
              clon(ifc),clat(ifc),fixvals(ifc))
    enddo
    
    
  !--------------------------------------------------------------
  !-- Print out a table listing of the locations of the minmax --
  !--------------------------------------------------------------
    if ( verb .ge. 3 ) then
      print *, ' '
      print *, '--------------------------------------------------'
      print *, " Fixes have been made for the object which"
      print *, " initial guess point is -lon: ", geslon, &
                  " -lat: ", geslat
      print *, " The centers are: "

      do ifc = 1, numfield%fixcenter
        print *, " For ", trim(fname%fixcenter(ifc)), "-", &
                    fname%fixcenter_lev(ifc)
        print *, " lon: ",glon(ptc_i(ifc))," lat: ",glat(ptc_(ifc)), &
                 " val: ", fixvals(ifc)
        endif
      enddo
      
      print *, " If values are -9999 it means a minmax wasn't found"
      
    endif


  !--------------------------------
  !-- Calculating the new center --
  !--------------------------------
    ! First, a mean of all the centers are done
    ! Then, it's recalculated by giving more weight to centers closer 
    !   to the previously calculated mean position
    ! If centers appeared to be at more than rad of the first mean
    !   position, they are not considered for the second one
    ! Note: we take a simple lon/lat average, without considering
    !   the Earth curvature
    ! A center will always be computed, even if it's the one provided


    ! simple average
    fix_lon = geslon
    fix_lat = geslat
    
    nfc = 0
    do ifc = 1, numfield%fixcenter
      if (fixvals(ifc) > -9998.) then
        nfc = nfc + 1
        fix_lat = fix_lat + glat(ptc_(ifc))
        if(fileinfo%grid_type == 'global') then
          if(glon(ptc_i(ifc)) - fix_lon > 180. ) then ! wrapping
            fix_lon = fix_lon + glon(ptc_i(ifc)) - 360
          else if (glon(ptc_i(ifc)) - fix_lon < -180. ) then ! wrapping
            fix_lon = fix_lon + glon(ptc_i(ifc)) + 360
          else
            fix_lon = fix_lon + glon(ptc_i(ifc))
          endif
        else
          fix_lon = fix_lon + glon(ptc_i(ifc))
        endif
      endif
    enddo
    
    fix_lon1 = fix_lon/(1+nfc)
    fix_lat1 = fix_lat/(1+nfc)


    ! weighted average
    call calcdist(fix_lon1, fix_lat1, geslon, geslat, &
            dist, degrees)
              
    if (dist > rad) then ! priority to the detection field
      fix_lon1 = geslon
      fix_lat1 = geslat
      wt = 1
    else
      wt = ((dist-rad)/rad)**2
    endif
      
    wts = wt
    fix_lon = geslon * wt
    fix_lat = geslat * wt
      
    
    nfc = 0
    do ifc = 1, numfield%fixcenter
      if (fixvals(ifc) > -9998.) then
      
        call calcdist(fix_lon1, fix_lat1, glon(ptc_i(ifc)), &
                glat(ptc_(ifc)), dist, degrees)
        
        if (dist > rad) then ! the parameter is not considered
          if (verb >= 3) then
            print *, " The variable ", trim(fname%fixcenter), "-", &
                    fname%fixcenter_lev
            print *, " appears to be too far from the others: ", dist
            print *, " it is discarded"
          endif
          fixvals(ifc) = -9999.
          clon(ifc) = -9999.
          clat(ifc) = -9999.
          
          exit
        endif
        
        wt = ((dist-rad)/rad)**2
        wts = wts + wt
        
        fix_lat = fix_lat + glat(ptc_(ifc))*wt
        
        if(fileinfo%grid_type == 'global') then
          if(glon(ptc_i(ifc)) - fix_lon > 180. ) then ! wrapping
            fix_lon = fix_lon + (glon(ptc_i(ifc)) - 360)*wt
          else if (glon(ptc_i(ifc)) - fix_lon < -180. ) then ! wrapping
            fix_lon = fix_lon + (glon(ptc_i(ifc)) + 360)*wt
          else
            fix_lon = fix_lon + glon(ptc_i(ifc))
          endif
        else
          fix_lon = fix_lon + glon(ptc_i(ifc))
        endif
        
      endif
    enddo
    
    fix_lon = fix_lon/wts
    fix_lat = fix_lat/wts

    
    
  !---------------------------------------------------------------------
  !-- Finally, determine the points to be masked, as being under the ---
  !---- influence of the object ----------------------------------------
  !---------------------------------------------------------------------
    ! if all the centers are in mask_object_min, mask_object is set to
    !   that one
    ! else if all the centers are in mask_object_max, the contour max
    !   (not closed) is the maximum value of the centers (in field detection)
    ! else, we need to relaunch check_contour with a contour_min fixed 
    !   to the maximum value of the centers, focing the calculation of
    !   the masks
    
    contour_max = fxy(pt_i, pt_j)
    in_mask = .true.
    
    do ifc = 1,numfield%fixcenter
      if (.not. mask_object_min(ptc_i(ifc), ptc_j(ifc))) then
        val = data_detection(ptc_i(ifc), ptc_j(ifc))
        
        if (fname%detec_minmax = 'min') then
          contour_max = max(contour_max, val)
        else
          contour_max = min(contour_max, val)
        endif
        
        if(.not. mask_object_max(ptc_i(ifc), ptc_j(ifc))) then
          in_mask = .false.
        endif
        
      endif
    enddo
    
    if(.not. in_mask) then
    
      call check_contour(pt_i, pt_j, data_detection, contour_max, &
              .true., fname%detec_minmax, mask_object_notused, &
              mask_object, contour_max, rcontour_max,cccret)
              
    else if (contour_max /= fxy(pt_i, pt_j)) then
    
      mask_object = mask_object_max
      do i=1,imax
      do j=1,jmax
        if(mask_object(i,j) .and. fxy(i,j) > contour_max) then
          mask_object(i,j) = .false.
        endif
      enddo
      enddo
      
      !we need to recompute rcontour_max
      nb = 0
      rcontour_max = 0
      do i = 1,imax
        do j = 1,jmax
          if(mask_object(i,j)) then
          
            call get_neighbours(i, j, fileinfo%grid_type, &
                imax, jmax, ptn, ptn_i, ptn_j)
            
            on_border = .false.
            do iptn=1,ptn
              if(.not. mask_object(ptn_i(iptn),ptn_j(iptn))) then
                on_border = .true.
              endif
            enddo
            
            if(on_border .and. (contour_max - fxy(i,j))*mm < 0) then
              call calcdist(glon(ptx_i), glat(ptx_j), glon(i), glat(j), &
                      dist,degrees)
              nb = nb + 1
              rcontour_max = rcontour_max + dist
            endif
            
          endif
        enddo
      enddo

      rcontour_max = rcontour_max/nb
      
    endif

 
  
  end subroutine
  
  
               
  
  
  
  !---------------------------------------------------------------------------------
  !----- Subroutine searching for new local min / max in the "detection" field -----
  !---------------------------------------------------------------------------------
  subroutine first_ges_center (ilast_object)
  !(its,masked_out,stormct,contour_info,maxmini,maxminj,ifgcret)
  
  ! imax,jmax,dx,dy in module data_param
  ! cparm, fxy -> it's data_detection
  ! cmaxmin -> it's fname%detec_minmax
  ! trkrinfo -> in module namelist
  ! ilast_object previously stormct
  
  
    implicit none
    
    integer ibeg, iend, jbeg, jend
    logical mask_detec(imax,jmax), mask_object(imax,jmax)
    integer point(2)
    real plastbar,rlastbar
    integer cccret
    
    
    
    if ( verb .ge. 2 ) then
      print *,' '
      print *,'*----------------------------*'
      print *,'* At top of first_ges_center *'
      print *,'*----------------------------*'

    endif
    
    
    
  !-------------------------------------------------    
  !-- Determine the indices of the box boundaries --
  !-------------------------------------------------
    if (trkrinfo%northbd < -998.0 .or. trkrinfo%southbd < -998.0 .or. &
         trkrinfo%westbd < -998.0  .or. trkrinfo%eastbd < -998.0) then
      ! The user did not specify a subgrid, so scan the whole domain
      ibeg = 1
      iend = imax
      jbeg = 1
      jend = jmax
    else    
  
      jbeg = max(1, int(((glatmax + dy - trkrinfo%northbd) / dy)))
      jend = min(jmax, int(((glatmax + 2*dy - trkrinfo%southbd) / dy)))

      ! To deal with boxes across the longitude cut (at 0°, 180°, or other)
      if (fileinfo%grid_type == 'global' .and. &
          trkrinfo%westbd .lt. glonmin) then
        ibeg = int(((trkrinfo%westbd - glonmin)  / dx) )
      else
        ibeg = max(1, int(((trkrinfo%westbd - glonmin + dx) / dx) ))
      endif

      if (fileinfo%grid_type == 'global' .and. &
          trkrinfo%eastbd .gt. glonmax) then
        iend = int(((trkrinfo%eastbd - glonmin + dx)  / dx) )
      else
        iend = min(imax, int(((trkrinfo%eastbd - glonmin + 2*dx)/dx)))
      endif
      
      if ((ibeg + imax) == iend) then !All the latitude are actually considered => simplification
        ibeg = 1
        iend = imax
      endif
      
    endif
    
    if ( verb .ge. 2 ) then
      print *, "The indies of the boundary box are :"    
      print *, "jbeg: ", jbeg
      print *, "jend: ", jend
      print *, "ibeg: ", ibeg
      print *, "iend: ", iend
    endif


  !-------------------------
  !-- define mask_detec --
  !------------------------- 
    ! hide the areas where the objects already being tracked are
    mask_detec = .not. masked_out

    ! hide the areas beyond the box of search
    if (ibeg < 1) then
      mask_detec((iend+1):(ibeg-1),) = .false.
    else if (iend > imax) then
      if (ibeg > imax + 1) then ! This case is only possible when grid is global and glonmin<0 
        mask_detec(1:(ibeg-1-imax),) = .false.
        mask_detec((iend+1-imax):imax,) = .false. ! we suppose glonmax>0 and thus iend<2*imax
      else
        mask_detec((iend+1-imax):(ibeg-1),) = .false.   
      endif
    else ! no wrapping case
      if (ibeg > 1) then
        mask_detec(1:(ibeg-1),)    = .false.
      endif
      if (iend < imax) then
        mask_detec((iend+1):imax,) = .false.
      endif
    endif
    
    if (jbeg > 1) then
      mask_detec(,1:(jbeg-1))    = .false.
    endif
    if (jend < jmax) then
      mask_detec(,(jend+1):jmax) = .false.
    endif

    ! hide the border of the domain
    mask_detec(,1)    = .false.
    mask_detec(,jmax) = .false.
    if (fileinfo%grid_type == 'regional') then
      mask_detec(1,)    = .false.
      mask_detec(imax,) = .false.
    endif

? old_ilast_object = ilast_object


  !---------------
  !-- Main loop --
  !---------------
    still_finding_valid_maxmins = .true.
    search_loop: do while (still_finding_valid_maxmins)
    
      if (fname%detec_minmax == 'min') then
        ptx = minloc(data_detection, mask = mask_detec)
        value = data_detection(point(1), point(2))
     
        if (value > trkrinfo%detec_thresh) then
          still_finding_valid_maxmins = .false.
          exit search_loop
        endif
                
      else
        ptx = maxloc(data_detection, mask = mask_detec)
        value = data_detection(point(1), point(2))
        
        if (value < trkrinfo%detec_thresh) then
          still_finding_valid_maxmins = .false.
          exit search_loop
        endif

      endif

!?      pt_ip1 = ptx(1) + 1
!?      pt_im1 = ptx(1) - 1
!?      pt_jp1 = ptx(2) + 1
!?      pt_jm1 = ptx(2) - 1
      ! wrapping, the definition of mask_detec made the cases only possible when the grid is global
!?      if (ptx(1) == 1)    pt_im1 = imax
!?      if (ptx(1) == imax) pt_ip1 = 1
    
      ! Check if it is a local minmax (the mask could have hidden one of the neighboring points)
!?      if (trkrinfo%detec_minmax == 'min') then
!?        if (value <= data_detection(pt_ip1,ptx(2)) .and. &
!?            value <= data_detection(ptx(1),pt_jm1) .and. &
!?            value <= data_detection(pt_im1,ptx(2)) .and. &
!?            value <= data_detection(ptx(1),pt_jp1)) then
!?          rough_gradient_check_okay = .true.
!?        else
!?          rough_gradient_check_okay = .false.
!?        endif
!?      else
!?        if (value >= data_detection(pt_ip1,ptx(2)) .and. &
!?            value >= data_detection(ptx(1),pt_jm1) .and. &
!?            value >= data_detection(pt_im1,ptx(2)) .and. &
!?            value >= data_detection(ptx(1),pt_jp1)) then
!?          rough_gradient_check_okay = .true.
!?        else
!?          rough_gradient_check_okay = .false.
!?        endif
!?      endif

!?      if (rough_gradient_check_okay) then
      if ( verb .ge. 3 ) then
        print *,'Found a possible max/min at ptx(1)= ', &
                  ptx(1),' ptx(2)= ',ptx(2)
      endif
      
    if (fname%detec_minmax == 'min') then
      contour = xcentval + trkrinfo%contint
    else
      contour = xcentval - trkrinfo%contint
    endif
          
      ! The point is considered if a closed contour is found around it
      call check_contour(ptx(1), ptx(2), data_detection, contour, &
              .false., fname%detec_minmax, mask_object_max, &
              mask_object_min, contour_max, rcontour_max,cccret)
              
              
      ! if no closed contour
      if (cccret /= 0) then
        if ( verb .ge. 3 ) then
          print *, "No contour has been found: no new object"
        endif        

        ! We don't want to find this local minimum or its 8 surrounding
        ! points again in a search on a subsequent iteration of this loop.
        call get_neighbours(ptx(1), ptx(2), fileinfo%grid_type, &
                imax, jmax, ptn, ptn_i, ptn_j)
        do iptn = 1,ptn
          mask_detec(ptn_i(iptn),ptn_j(iptn)) = .false.
        enddo
        mask_detec(ptx(1),ptx(2)) = .false.
          
        cycle search_loop

      endif 


      ! A closed contour has been found
      ! But there is no more space in the object arrays
      if (ilast_object >= max_object) then
        if ( verb .ge. 1 ) then
          print *, " "
          print *, "ERROR 51 in first_ges_center:"
          print *, "The maximum number of object has been reached"
          print *, "The program can be relauched from "
          print *, "that timestep; ts_tim(its): ", ts_tim(its)
          print *, "Using the record of the last objects being"
          print *, "tracked will discard all the old storms"
        endif  
        STOP 51
      endif        
          
      ilast_object = ilast_object + 1
      if ( verb .ge. 3 ) then
        print *, "A contour has been found: this is a new object"
      endif
      
      ! name the object
      call fixcenter (ptx(1), ptx(2), mask_object_max, &
                mask_object_min, mask_object, fixlon, fixlat, fixvals, &
                contour_max, rcontour_max)

      call output...

???? !!!!! Need to fill the object record here using ptx(1) and ptx(2) (maxmini, maxminj)
? ! do the fptx(1)ing of the object center here, so that mask_detec (masked_out) can be updated
?
?

      ! fill mask_detec with new false and masked_out with true
      mask_detec = mask_detec .and. .not. mask_object_max
      masked_out   = masked_out .or. mask_object

    enddo search_loop

! Print the number of new lows found, rise a warning if none




! filling the object record --> could be done from the implementation of find_all_maxmins above
    if (old_last_storm > last_storm)
        do n = isstart,stormct
          if (trkrinfo%type == 'midlat') then
            storm(n)%tcv_center = 'MIDL'
          else if (trkrinfo%type == 'tcgen') then
            storm(n)%tcv_center = 'TCG '
          endif
          slonfg(n,ifh) = glonmin + (maxmini(n)-1)*dx
          slatfg(n,ifh) = glatmax - (maxminj(n)-1)*dy
          storm(n)%tcv_stspd = -99
          storm(n)%tcv_stdir = -99
          write (storm(n)%tcv_storm_id,'(i4.4)') n
          write (storm(n)%tcv_storm_name,'(i4.4)') n
          stormswitch(n) = 1
          if (cparm == 'mslp') then

            if ( verb .ge. 3 ) then
              write (6,71) maxmini(n),maxminj(n),slonfg(n,ifh)
     &             ,360.-slonfg(n,ifh),slatfg(n,ifh)
     &             ,slp(maxmini(n),maxminj(n))/100.0
            endif

          endif
        enddo  
  
  
  
  
  end subroutine
  

end module




