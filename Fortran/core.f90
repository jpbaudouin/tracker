!#######################################################################
!########## DEFINITION OF THE CORE PARAMETERS AND SUBROUTINES ##########
!#######################################################################


!25-10-2017 Creation
!12-12-2017 1st Commit -> pi, dtr
!14-12-2017 2nd commit -> check_closed_contour, first_ges_center

! Is there a better way to change max_object than recompile it?
?? ! Watch out: mask_detec and masked_out have opposite meaning for their values
?? ! all .true. are considered for mask_detec, while all .false. are for masked_out


! CONTAINS :
!   The modules core with different parameters




module core
  ! core subroutine and parameters
  
  use data_param
  use namelist_trk

  implicit none
    
  real, parameter :: pi = 4. * atan(1.)
  real, parameter :: dtr = 4. * atan(1.)/180.0 !pi/180

  integer, parameter :: max_object = 5000
  
?!  mask_out and ilast_object
  
  
  contains

  !---------------------------------------------------------------------------------
  !----------- Subroutine searching for a closed contour around a point ------------
  !---------------------------------------------------------------------------------

  subroutine check_closed_contour (ptx_i,ptx_j,fxy,cmaxmin, &
               multi_contour,mask_object,plastbar,rlastbar,ret)

! masked_out might be an argument if different than "area influenced by other object"
! It seems it's values are updated in the subroutine'
! merge closed_contour with ret (icccret)
! num_requested_contours modified into logical multi_contour (not implemented yet)
! get_last_isobar_flag is supressed


    implicit none
    
    integer, intent(in) :: ptx_i, ptx_j
    real,    intent(in) :: fxy(:,:)
    logical, intent(in) :: multi_contour
    character*3, intent(in) :: cmaxmin
    
    logical, intent(out) :: mask_object(:,:)    
    real,    intent(out) :: plastbar, rlastbar
    integer, intent(out) :: ret = 0

    real :: xcentval = fxy(ix,jx)
    real contour, dist, degree
    integer ipt, pt_max, iptb, iptn
    integer pt_ip1, pt_im1, pt_jp1, pt_jm1
    integer pt_i(imax*jmax), pt_j(imax*jmax), ptn_i(8), ptn_j(8)
    logical contour_close
ptn_i
      
    if ( verb .ge. 3 ) then
      print *,' '
      print *,'*-----------------------------------------------*'
      print *,'* Top of check_closed_contour, ix= ',ix,' jx= ',jx
      print *,'*-----------------------------------------------*'
      print *,' '
      print *,'fxy(ix,jx)= ',fxy(ix,jx)
    endif


    ! set up the contour interval that will be used, from the value of the central point 
    if (cmaxmin == 'min') then
      contour = xcentval + trkrinfo%contint
    else
      contour = xcentval - trkrinfo%contint
    endif

! In case of multi_contour = 'y' this loop need to be implemented
!      successive_contours_loop: do while (num_found_contours <
!     &     num_requested_contours)


    !initialise the search to the first point
    ipt = 1 ! current index
    pt_max = 1 ! maximum index
    pt_i(1) = ptx_i ! the points to be treated
    pt_j(1) = ptx_j 
    mask_object(ptx_i, ptx_j) = .false.
    
    iptb = 0 ! the number of points at the border of the contour
          
    ! Loop on the points (pt_i, pt_j) 
    scan_contour: do while (ipt <= pt_max)
    
      ! if we hit the border of a regional grid, there is no contour
      if (fileinfo%grid_type == 'regional' .and. &
         (pt_i(ipt) == 1 .or. pt_i(ipt) == imax .or. &
          pt_j(ipt) == 1 .or. pt_j(ipt) == jmax)) then
        if ( verb .ge. 3 ) then
          print *,'We reach the border of a regional grid'
          print *,'There is no close contour found'
        endif
        contour_close = .false.
          ret = 61
        exist scan_contour
      endif          

      ! Determine i / j of the 8 neighboring points
      pt_ip1 = pt_i(ipt) + 1
      pt_im1 = pt_i(ipt) - 1
      pt_jp1 = pt_j(ipt) + 1
      pt_jm1 = pt_j(ipt) - 1
      if (fileinfo%grid_type == 'global') then ! wrapping
        if (pt_i(ipt) == 1)    pt_im1 = imax
        if (pt_i(ipt) == imax) pt_ip1 = 1
      endif

      ! list of the 8 neighboring points
      ptn_i = [pt_ip1, pt_ip1   , pt_ip1, pt_i(ipt), &
               pt_im1, pt_im1   , pt_im1, pt_i(ipt)  ]
      ptn_j = [pt_jp1, pt_j(ipt), pt_jm1, pt_jm1,    &
               pt_jm1, pt_j(ipt), pt_jp1, pt_jp1     ]
      ptn = 8

      ! if we hit the Nth/ Sth border of a global grid, we continue the search 
      !  for contour with the other neighboring points
      if (fileinfo%grid_type == 'global') then
        if (pt_j(ipt) == 1) then
          ptn_i = [pt_ip1, pt_ip1   , pt_im1   , pt_im1, pt_i(ipt)]
          ptn_j = [pt_jp1, pt_j(ipt), pt_j(ipt), pt_jp1, pt_jp1   ]
          ptn = 5
        else if (pt_j(ipt) == jmax) then
          ptn_i = [pt_ip1   , pt_ip1, pt_i(ipt), pt_im1, pt_im1   ]
          ptn_j = [pt_j(ipt), pt_jm1, pt_jm1   , pt_jm1, pt_j(ipt)]
          ptn = 5
      endif
      
      scan_neighbour : do iptn=1,ptn
        if (masked_out(ptn_i(iptn), ptn_j(iptn))) then
          if ( verb .ge. 3 ) then
            print *,'We reach the border of another object'
            print *,'There is no close contour found'
          endif
          contour_close = .false.
          ret = 62
          exist scan_contour ! It seems to be possible to exit two do loop like this

          
        else if (cmaxmin == 'min' .and. &
                 fxy(ptn_i(iptn), ptn_j(iptn)) < xcentval) then
          if ( verb .ge. 3 ) then
            print *, "We reach a point with a lower value than ", &
                      "the limit"
            print *, "fxy(i,j) :" fxy(ptn_i(iptn), ptn_j(iptn))
            print *, "limit: ", xcentval
            print *, "There is no close contour found"
          endif
          contour_close = .false.
          ret = 63
          exist scan_contour

        else if (cmaxmin == 'max' .and. &
                 fxy(ptn_i(iptn), ptn_j(iptn)) > xcentval) then
          if ( verb .ge. 3 ) then
            print *, "We reach a point with a higher value than ", &
                      "the limit"
            print *, "fxy(i,j) :" fxy(ptn_i(iptn), ptn_j(iptn))
            print *, "limit: ", xcentval
            print *, "There is no close contour found"
          endif
          contour_close = .false.
          ret = 64
          exist scan_contour

          
        else if (.not. mask_object(ptn_i(iptn), ptn_j(iptn))) then
          if ((cmaxmin == 'max' .and. &
                fxy(ptn_i(iptn), ptn_j(iptn)) >= contour) .or.
              (cmaxmin == 'min' .and. &
                fxy(ptn_i(iptn), ptn_j(iptn)) <= contour)) then
        
            mask_object(ptn_i(iptn), ptn_j(iptn)) = .true.
            pt_max = pt_max + 1
            pt_i(pt_max) = ptn_i(iptn)
            pt_j(pt_max) = ptn_j(iptn)
            
          else ! the point from which derive the neighbours is on the border of the contour
            iptb = iptb + 1
            call calcdist (glon(ptx_i),glat(ptx_j), &
                   glon(pt_i(ipt)),glat(pt_j(ipt)),dist,degrees)
            rlastbar = rlastbar + dist
          endif
          
        endif
      enddo scan_neighbour

      ipt = ipt + 1
      
    enddo scan_contour

    if (contour_close) then
      plastbar = contour
      rlastbar = rlastbar / float(iptb)
    else
      plastbar = -999.
      rlastbar = -999.
    endif

  end subroutine
  
  
  
  
  !---------------------------------------------------------------------------------
  !----- Subroutine searching for new local min / max in the "detection" field -----
  !---------------------------------------------------------------------------------
  subroutine first_ges_center (ilast_object)
  !(its,masked_out,stormct,contour_info,maxmini,maxminj,ifgcret)
  
  ! imax,jmax,dx,dy in module data_param
  ! cparm, fxy -> it's data_detection
  ! cmaxmin -> it's trkrinfo%detec_minmax
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
  !-- Determining the indices of the box boundary --
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


? old_ilast_object = ilast_object




  !-------------------------
  !-- defining mask_detec --
  !------------------------- 
    ! hide the areas where the objects already being tracked are
    mask_detec = .not. masked_out

    ! hide the areas beyond the box of search
    if (ibeg .le. 0) then
      mask_detec((iend+1):(ibeg-1),) = .false.
    else if (iend .gt. imax) then
      if (ibeg .gt. imax + 1) then ! This case is only possible when grid is global and glonmin<0 
        mask_detec(1:(ibeg-1-imax),) = .false.
        mask_detec((iend+1-imax):imax,) = .false. ! we suppose glonmax>0 and thus iend<2*imax
      else
        mask_detec((iend+1-imax):(ibeg-1),) = .false.   
      endif
    else ! no wrapping case
      if (ibeg .gt. 1) then
        mask_detec(1:(ibeg-1),)    = .false.
      endif
      if (iend .lt. iend) then
        mask_detec((iend+1):imax,) = .false.
      endif
    endif
    
    if (jbeg .gt. 1) then
      mask_detec(,1:(jbeg-1))    = .false.
    endif
    if (jend .lt. jmax) then
      mask_detec(,(jend+1):jmax) = .false.
    endif

    ! hide the border of the domain
    mask_detec(,1)    = .false.
    mask_detec(,jmax) = .false.
    if (fileinfo%grid_type == 'regional') then
      mask_detec(1,)    = .false.
      mask_detec(imax,) = .false.
    endif


  !---------------
  !-- Main loop --
  !---------------
    still_finding_valid_maxmins = .true.
    search_loop: do while (still_finding_valid_maxmins)
    
      if (trkrinfo%detec_minmax == 'min') then
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
      
??????! Is This part so useful ?? check_closed_contour will find out if this not a local maxmin
?      ! The four neighboring points
?      pt_ip1 = ptx(1) + 1
?      pt_im1 = ptx(1) - 1
?      pt_jp1 = ptx(2) + 1
?      pt_jm1 = ptx(2) - 1
?      ! wrapping, the definition of mask_detec made the cases only possible when the grid is global
?      if (ptx(1) == 1)    pt_im1 = imax
?      if (ptx(1) == imax) pt_ip1 = 1
?        
?      ! Check if it is a local minmax (the mask could have hidden one of the neighboring points)
?      if (trkrinfo%detec_minmax == 'min') then
?        if (value <= data_detection(pt_ip1,ptx(2)) .and. &
?            value <= data_detection(ptx(1),pt_jm1) .and. &
?            value <= data_detection(pt_im1,ptx(2)) .and. &
?            value <= data_detection(ptx(1),pt_jp1)) then
?          rough_gradient_check_okay = .true.
?        else
?          rough_gradient_check_okay = .false.
?        endif
?      else
?        if (value >= data_detection(pt_ip1,ptx(2)) .and. &
?            value >= data_detection(ptx(1),pt_jm1) .and. &
?            value >= data_detection(pt_im1,ptx(2)) .and. &
?            value >= data_detection(ptx(1),pt_jp1)) then
?          rough_gradient_check_okay = .true.
?        else
?          rough_gradient_check_okay = .false.
?        endif
?      endif
?
?      if (rough_gradient_check_okay) then
        if ( verb .ge. 3 ) then
          print *,'Found a possible max/min at ptx(1)= ', &
                    ptx(1),' ptx(2)= ',ptx(2)
        endif


        ! The point is considered if a closed contour is found around it
        call check_closed_contour(ptx(1), ptx(2), data_detection, &
               masked_out, trkrinfo%detec_minmax, .false., &
               mask_object,plastbar,rlastbar,cccret)

          
        ! if no closed contour
        if (cccret /= 0) then
          if ( verb .ge. 3 ) then
            print *, "No contour has been found: no new object"
          endif        

          ! We don't want to find this local minimum or its 8 surrounding
          ! points again in a search on a subsequent iteration of this loop.
          mask_detec(ptx(1),ptx(2)) = .false.
          mask_detec(ptx(1),pt_jp1) = .false.
          mask_detec(pt_ip1,pt_jp1) = .false.
          mask_detec(pt_ip1,ptx(2)) = .false.
          mask_detec(pt_ip1,pt_jm1) = .false.
          mask_detec(ptx(1),pt_jm1) = .false.
          mask_detec(pt_im1,pt_jm1) = .false.
          mask_detec(pt_im1,ptx(2)) = .false.
          mask_detec(pt_im1,pt_jp1) = .false.   
          
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

        call fixcenter
        call check_closed_contour
        call output...

???? !!!!! Need to fill the object record here using ptx(1) and ptx(2) (maxmini, maxminj)
? ! do the fptx(1)ing of the object center here, so that mask_detec (masked_out) can be updated
?
?

        ! fill mask_detec with new false and masked_out with true


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




