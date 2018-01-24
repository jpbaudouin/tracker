!#######################################################################
!####### DEFINITION OF THE SUBROUTINES FOR DISTANCE CALCULATION ########
!#######################################################################

!17-01-2018 1st Commit -> creation
!24-01-2018 2nd commit -> working version



! CONTAINS :
!   The module "grid_projection". 
!   It includes some parameters for distance calcultation
!   It includes the subroutines:
!       - calc_dist_dir(lon1,lat1,lon2,lat2,calcdir,dist,dir)
!       - calc_next_point(lon1,lat1,dist,dir,lon2,lat2)




!###########################################################################################
!############## Module including subrtoutines related to the grid projection ###############
!###########################################################################################



module grid_projection


  implicit none

  real, parameter :: dtr = 4. * atan(1.)/180.0 !pi/180
  real, parameter :: erad = 6371.0 ! Earth's radius (km)


  contains
  

  !---------------------------------------------------------------------------------
  !------ Subroutine calculating distance and direction between two points ---------
  !---------------------------------------------------------------------------------
  subroutine calc_dist_dir(lon1,lat1,lon2,lat2,calcdir,dist,dir)

  ! The algorithm is based on a projection of the path over a spherical
  !    Earth, using the great-circle. It computes the distance between 
  !    two points and the direction at the last one (final bearing)
  ! The distance uses the haversine formula. Both distance and direction
  !    formula are adapted from:
  !    https://www.movable-type.co.uk/scripts/latlong.html

  
    ! INPUT
    ! lon1    longitude of the starting point
    ! lat1    latitude  of the starting point
    ! lon1    longitude of the ending point
    ! lat1    latitude  of the ending point
    ! calcdir indicate whether to compute the direction. It is only 
    !           needed in get_speed_dir for the extrapolation method

   
    ! OUTPUT
    ! dist    distance between the point (in km)
    ! dir     direction (in degree) in which the object is moving at the
    !           ending point


    use namelist_trk

    implicit none

    real, intent(in) :: lon1, lon2, lat1, lat2
    logical, intent(in) :: calcdir
    real, intent(out) :: dist, dir
      
    real phi1, phi2, lambda1, lambda2, a      
    
    if (lat1 > 90 .or. lat1 < -90 .or. lat2>90 .or. lat2>90) then
      if (verb >= 1) then
        print *, ""
        print *, "!!! ERROR 301: in calc_dist_dir"
        print *, "!!! One of the latitude exceeded 90 or -90°:"
        print *, "!!! lat1: ", lat1
        print *, "!!! lat2: ", lat2
      endif 
      STOP 301
    endif  
      
      
    phi1 = lat1*dtr
    phi2 = lat2*dtr
      
    lambda1 = lon1*dtr
    lambda2 = lon2*dtr
      

    a = cos(phi1)*cos(phi2)*sin((lambda2 - lambda1)/2)**2 + &
          sin((phi2 - phi1)/2)**2 
    dist = erad*2*atan2(sqrt(a), sqrt(1-a))
      
    if (calcdir) then
      ! This is the direction from where the object arrive, in radian
      dir = atan2(sin(lambda1 - lambda2)*cos(phi1), &
                  cos(phi2)*sin(phi1) - &
                    sin(phi2)*cos(phi1)*cos(lambda1-lambda2))
        
      dir = dir/dtr + 180 ! right direction in degree
      
      ! to keep the values between 0 and 360
      if (dir < 0) then
        dir = dir + 360
      endif
    
    else
      dir = -9999
    endif
      
         
      
  end subroutine
    


  !---------------------------------------------------------------------------------
  !----- Subroutine calculating the next point from distance and direction  --------
  !---------------------------------------------------------------------------------
  subroutine calc_next_point(lon1,lat1,dist,dir,lon2,lat2)

  ! Same argument as calc_dist_dir, but find the next point acording to
  !    a direction and a distance from a starting point.
  ! The formula are also adapted from the same website
  
  
    ! INPUT
    ! lon1    longitude of the starting point
    ! lat1    latitude  of the starting point
    ! dist    distance from the starting point of the ending point 
    !           we look for (in m)
    ! dir     direction (in degree) in which the object is moving at the
    !           starting point
    
    ! OUTPUT
    ! lon1    longitude of the ending point
    ! lat1    latitude  of the ending point
  
  
    implicit none 
    
    real, intent(in) :: lon1, lat1, dist, dir
    real, intent(out) :: lon2, lat2
    
    real delta, theta, phi1, lambda1, phi2, lambda2
  
    
    delta = dist/erad ! angular distance  
    theta = dir*dtr
    phi1  = lat1*dtr
    lambda1 = lon1*dtr

  
    phi2 = asin(sin(phi1)*cos(delta) + cos(phi1)*sin(delta)*cos(theta))
    lambda2 = lambda1 + atan2(sin(theta)*sin(delta)*cos(phi1), &
                              cos(delta) - sin(phi1)*sin(phi2))  
                            
    lat2 = phi2 / dtr
    lon2 = lambda2 / dtr
    
    if (lon2 < 0) then
      lon2 = lon2 + 360
    endif

    
    if (lon2 > 360) then
      lon2 = lon2 - 360
    endif
  
  
  end subroutine
  
    

end module





