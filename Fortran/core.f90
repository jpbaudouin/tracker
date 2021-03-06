!#######################################################################
!########## DEFINITION OF THE CORE PARAMETERS AND SUBROUTINES ##########
!#######################################################################


!25-10-2017 Creation
!12-12-2017 1st Commit -> pi, dtr
!14-12-2017 2nd commit -> check_contour, first_ges_center
!21-12-2017 3rd commit -> module no_dependance, get_neighbours
!                         get_ij_bounds, find_minmax, fixcenter
!                         correction check_contour (check_contour), first_ges_center
!12-01-2018 4th commit -> correction in get_neighbours
!                         writing check_contour, find_minmax and fixcenter help
!                         finishing first_ges_center
!15-01-2018 5th commit -> adding get_speed_dir and adding speed/direction in output (first_ges_center)
!24-01-2018 6th commit -> working version



! Is there a better way to change max_object than recompile it?

! Note: mask_detec and masked_out have opposite meaning for their values
! all .true. are considered for mask_detec, while all .false. are for masked_out


! CONTAINS :
!   The module "no_dependance" which includes subroutines without dependance :
!       - get_neighbours(pt_i, pt_j, grid, imax, jmax, ptn, ptn_i, ptn_j)
!       - qsort(x,ind,n)

!   The module "core_tracker" with few parameters and the subroutines:
!       - get_ij_bounds (geslon, geslat, rad, ibeg, jbeg, iend, jend, ret)
!       - check_contour (ptx_i, ptx_j, fxy, cminmax, contour_min, 
!                        forced, mask_object, mask_object_cc, 
!                        contour_max, rcontour_max,ret)
!       - find_minmax(geslon, geslat, fxy, minmax, ptc_i, ptc_j, 
!                     clon, clat, cval, found)
!       - fixcenter (pt_i, pt_j, mask_object_max, mask_object_min,
!                    mask_object, fix_lon, fix_lat, fixvals,
!                    contour_max, rcontour_max)
!       - get_speed_dir(iobj, its, fixlon, fixlat, mask_object, first_detection)  
!       - first_ges_center (its, ilast_object)




!###########################################################################################
!#################### Module including subrtoutines without dependances ####################
!###########################################################################################
!-----------------------------------------------------------------------------------

module no_dependance

  contains
  

  !---------------------------------------------------------------------------------
  !-------------- Subroutine returning a list of neighboring points ----------------
  !---------------------------------------------------------------------------------
  subroutine get_neighbours(pt_i, pt_j, grid, imax, jmax, ptn, ptn_i, ptn_j)
  
! This function return a list of ptn points represented by their 
!   indices (ptn_i, ptn_j) that are around a point of index (pt_i, pt_j).
! It considers a maximum of 8 neighboring points of (pt_i, pt_j).
! If the grid is "global", it considers a possble wrapping along i
! If the grid is "regional" and (pt_i, pt_j) is on its edge,
!   it will reduce ptn to the number of existing neighbours (5 or 3)
! There is no dependance on the original grid, in order to be used on 
!   local grid as given by get_ij_bounds

! previously get_ijplus1_check_wrap
  
  
  ! INPUT
  ! pt_i      Index along i of the central point (0 < pt_i <= i_max)
  ! pt_j      Index along j of the central point (0 < pt_j <= j_max)
  ! grid      the type of grid for possible wrapping, either "regional"
  !             or "global"
  ! imax      The maximum number of point along i on the grid (> 1)
  ! jmax      The maximum number of point along j on the grid (> 1)

  ! OUTPUT
  ! ptn       number of neighboring points (8, 5 or 3)
  ! ptn_i     List of index along i of the neighboring points
  ! ptn_j     List of index along j of the neighboring points
  
  
    implicit none
  
    integer, intent(in)  :: pt_i, pt_j, imax, jmax
    character*8, intent(in) :: grid
    integer, intent(out) :: ptn, ptn_i(8), ptn_j(8)
  
    integer pt_ip1, pt_im1, pt_jp1, pt_jm1
    
    
    ! error on the grid sent
    if (imax < 1 .or. jmax < 1) then
!      if ( verb .ge. 1 ) then
        print *, ""
        print *, "!!! ERROR 91 in get_neighbours:"
        print *, "!!! The grid provided is too small"
        print *, "!!! imax: ", imax, " jmax: ", jmax
!      endif
      STOP 91
    endif


  !-------------------------------------------------  
  !-- Determine i / j of the 8 neighboring points --
  !-------------------------------------------------
    pt_ip1 = pt_i + 1
    pt_im1 = pt_i - 1
    pt_jp1 = pt_j + 1
    pt_jm1 = pt_j - 1
    if (grid == 'global') then ! wrapping
      if (pt_i == 1)    pt_im1 = imax
      if (pt_i == imax) pt_ip1 = 1
    endif

    ! list of the 8 neighboring points
    ptn_i = [pt_ip1, pt_ip1, pt_ip1, pt_i, &
             pt_im1, pt_im1, pt_im1, pt_i]
    ptn_j = [pt_jp1, pt_j, pt_jm1, pt_jm1, &
             pt_jm1, pt_j, pt_jp1, pt_jp1]
    ptn = 8

  !---------------------------------------------------------------------
  !-- if the border of a grid is hit, only consider the reamining points
  !-- (if model is global, only Sth and Nth or considered) -------------
  !---------------------------------------------------------------------
    ! along j
    if (pt_j == 1) then
      ptn = 5
      ptn_i(1:ptn) = [pt_ip1, pt_ip1, pt_im1, pt_im1, pt_i]
      ptn_j(1:ptn) = [pt_jp1, pt_j,   pt_j,   pt_jp1, pt_jp1]
    else if (pt_j == jmax) then
      ptn = 5
      ptn_i(1:ptn) = [pt_ip1, pt_ip1, pt_i,   pt_im1, pt_im1]
      ptn_j(1:ptn) = [pt_j,   pt_jm1, pt_jm1, pt_jm1, pt_j]
    endif

    ! along i and corners
    if (grid == 'regional') then
      if (pt_i == 1) then
        if (pt_j == 1) then
          ptn = 3
          ptn_i(1:ptn) = [1,2,2]
          ptn_j(1:ptn) = [2,2,1]
        else if (pt_j == jmax) then
          ptn = 3
          ptn_i(1:ptn) = [1,2,2]
          ptn_j(1:ptn) = [jmax-1,jmax-1,jmax]
        else
          ptn = 5
          ptn_i(1:ptn) = [pt_ip1, pt_ip1, pt_ip1, pt_i,   pt_i]
          ptn_j(1:ptn) = [pt_jp1, pt_j,   pt_jm1, pt_jm1, pt_jp1]
        endif 
      else if (pt_i == imax) then
        if (pt_j == 1) then
          ptn = 3
          ptn_i(1:ptn) = [imax,imax-1,imax-1]
          ptn_j(1:ptn) = [2,2,1]
        else if (pt_j == jmax) then
          ptn = 3
          ptn_i(1:ptn) = [imax,imax-1,imax-1]
          ptn_j(1:ptn) = [jmax-1,jmax-1,jmax]
        else
          ptn = 5
          ptn_i(1:ptn) = [pt_i,   pt_im1, pt_im1, pt_im1, pt_i]
          ptn_j(1:ptn) = [pt_jm1, pt_jm1, pt_j,   pt_jp1, pt_jp1]
        endif
      endif
    endif
             
  end subroutine


  !---------------------------------------------------------------------------------
  !-------------------- Subroutine sorting values in a vector ----------------------
  !---------------------------------------------------------------------------------
  SUBROUTINE qsort(x,ind,n)

  ! Code converted using TO_F90 by Alan Miller
  ! Date: 2002-12-18  Time: 11:55:47
  ! ********************************************************************
  ! 
  !                                                         ROBERT RENKA
  !                                                 OAK RIDGE NATL. LAB.
  ! 
  ! THIS SUBROUTINE USES AN ORDER N*LOG(N) QUICK SORT TO SORT A REAL 
  !   (dp) ARRAY X INTO INCREASING ORDER.  THE ALGORITHM IS AS FOLLOWS.  
  !   IND IS INITIALIZED TO THE ORDERED SEQUENCE OF INDICES 1,...,N, AND 
  !   ALL INTERCHANGES ARE APPLIED TO IND.  X IS DIVIDED INTO TWO 
  !   PORTIONS BY PICKING A CENTRAL ELEMENT T.  THE FIRST AND LAST 
  !   ELEMENTS ARE COMPARED WITH T, AND INTERCHANGES ARE APPLIED AS 
  !   NECESSARY SO THAT THE THREE VALUES ARE IN ASCENDING ORDER.  
  !   INTERCHANGES ARE THEN APPLIED SO THAT ALL ELEMENTS GREATER THAN T
  !   ARE IN THE UPPER PORTION OF THE ARRAY AND ALL ELEMENTS LESS THAN T
  !   ARE IN THE LOWER PORTION.  THE UPPER AND LOWER INDICES OF ONE OF
  !   THE PORTIONS ARE SAVED IN LOCAL ARRAYS, AND THE PROCESS IS 
  !   REPEATED ITERATIVELY ON THE OTHER PORTION.  WHEN A PORTION IS
  !   COMPLETELY SORTED, THE PROCESS BEGINS AGAIN BY RETRIEVING THE 
  !   INDICES BOUNDING ANOTHER UNSORTED PORTION.
  

    ! INPUT 
    ! N     LENGTH OF THE ARRAY X.
    ! X     VECTOR OF LENGTH N TO BE SORTED.

    ! OUTPUT 
    ! IND   SEQUENCE OF INDICES 1,...,N PERMUTED IN THE SAME FASHION AS 
    !         X WOULD BE. THE ORDERING ON X IS DEFINED BY Y(I) = X(IND(I))

    ! NOTE: IU AND IL MUST BE DIMENSIONED >= LOG(N) WHERE LOG HAS BASE 2.
  
    
      IMPLICIT NONE
      INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

      REAL, INTENT(IN)  :: x(n)
      INTEGER, INTENT(OUT)   :: ind(n)
      INTEGER, INTENT(IN)    :: n




      INTEGER   :: iu(21), il(21)
      INTEGER   :: m, i, j, k, l, ij, it, itt, indx
      REAL      :: r
      REAL (dp) :: t

      ! LOCAL PARAMETERS -

      ! IU,IL =  TEMPORARY STORAGE FOR THE UPPER AND LOWER
      !            INDICES OF PORTIONS OF THE ARRAY X
      ! M =      INDEX FOR IU AND IL
      ! I,J =    LOWER AND UPPER INDICES OF A PORTION OF X
      ! K,L =    INDICES IN THE RANGE I,...,J
      ! IJ =     RANDOMLY CHOSEN INDEX BETWEEN I AND J
      ! IT,ITT = TEMPORARY STORAGE FOR INTERCHANGES IN IND
      ! INDX =   TEMPORARY INDEX FOR X
      ! R =      PSEUDO RANDOM NUMBER FOR GENERATING IJ
      ! T =      CENTRAL ELEMENT OF X

      IF (n <= 0) RETURN

      ! INITIALIZE IND, M, I, J, AND R

      DO  i = 1, n
        ind(i) = i
      END DO
      m = 1
      i = 1
      j = n
      r = .375

      ! TOP OF LOOP

   20 IF (i >= j) GO TO 70
      IF (r <= .5898437) THEN
        r = r + .0390625
      ELSE
        r = r - .21875
      END IF

      ! INITIALIZE K

   30 k = i

      ! SELECT A CENTRAL ELEMENT OF X AND SAVE IT IN T

      ij = i + r*(j-i)
      it = ind(ij)
      t = x(it)

      ! IF THE FIRST ELEMENT OF THE ARRAY IS GREATER THAN T,
      !   INTERCHANGE IT WITH T

      indx = ind(i)
      IF (x(indx) > t) THEN
        ind(ij) = indx
        ind(i) = it
        it = indx
        t = x(it)
      END IF

      ! INITIALIZE L

      l = j

      ! IF THE LAST ELEMENT OF THE ARRAY IS LESS THAN T,
      !   INTERCHANGE IT WITH T
      indx = ind(j)
      IF (x(indx) >= t) GO TO 50
      ind(ij) = indx
      ind(j) = it
      it = indx
      t = x(it)

      ! IF THE FIRST ELEMENT OF THE ARRAY IS GREATER THAN T,
      !   INTERCHANGE IT WITH T

      indx = ind(i)
      IF (x(indx) <= t) GO TO 50
      ind(ij) = indx
      ind(i) = it
      it = indx
      t = x(it)
      GO TO 50

      ! INTERCHANGE ELEMENTS K AND L

   40 itt = ind(l)
      ind(l) = ind(k)
      ind(k) = itt

      ! FIND AN ELEMENT IN THE UPPER PART OF THE ARRAY WHICH IS
      !   NOT LARGER THAN T

   50 l = l - 1
      indx = ind(l)
      IF (x(indx) > t) GO TO 50

      ! FIND AN ELEMENT IN THE LOWER PART OF THE ARRAY WHCIH IS NOT SMALLER THAN T

   60 k = k + 1
      indx = ind(k)
      IF (x(indx) < t) GO TO 60

      ! IF K <= L, INTERCHANGE ELEMENTS K AND L

      IF (k <= l) GO TO 40

      ! SAVE THE UPPER AND LOWER SUBSCRIPTS OF THE PORTION OF THE
      !   ARRAY YET TO BE SORTED

      IF (l-i > j-k) THEN
        il(m) = i
        iu(m) = l
        i = k
        m = m + 1
        GO TO 80
      END IF

      il(m) = k
      iu(m) = j
      j = l
      m = m + 1
      GO TO 80


      ! BEGIN AGAIN ON ANOTHER UNSORTED PORTION OF THE ARRAY

   70 m = m - 1
      IF (m == 0) RETURN
      i = il(m)
      j = iu(m)

   80 IF (j-i >= 11) GO TO 30
      IF (i == 1) GO TO 20
      i = i - 1

      ! SORT ELEMENTS I+1,...,J.  NOTE THAT 1 <= I < J AND J-I < 11.

   90 i = i + 1
      IF (i == j) GO TO 70
      indx = ind(i+1)
      t = x(indx)
      it = indx
      indx = ind(i)
      IF (x(indx) <= t) GO TO 90
      k = i

  100 ind(k+1) = ind(k)
      k = k - 1
      indx = ind(k)
      IF (t < x(indx)) GO TO 100

      ind(k+1) = it
      GO TO 90
      END SUBROUTINE qsort
  
  
  
end module


!###########################################################################################
!#################### Module including core subrtoutines and parameters ####################
!###########################################################################################
!----------------------------------------------------------------------------------

module core_tracker
  ! core subroutine and parameters
  
  use data_param
  use namelist_trk
  use write_read_output
  use grid_projection
  use no_dependance

  implicit none
  
  integer, parameter :: max_object = 9999
  logical, save, allocatable :: masked_out(:,:)
  logical beingTracked(max_object) ! previously stormswitch, true if a storm needs to be tracked
 
  
  contains
  


  !---------------------------------------------------------------------------------
  !-------------- Subroutine reducing array size for local analysis ----------------
  !---------------------------------------------------------------------------------

  subroutine get_ij_bounds (geslon, geslat, rad, &
               ibeg, jbeg, iend, jend, ret)
  
! It returns the indices of a subgrid that includes all the points
!   within a distance rad of the point of coordinate (geslon, geslat)
! Note: the point (geslon, geslat) doesn't need to be on the main grid
! Note: the points at the edge of this subgrid are all farther than
!     the distance rad of the central point.
!   Effectively, ibeg, jbeg, iend, jend represent the indeces of the 
!     first point farther than rad towards, respectively the West, 
!     the North, the East and the South of the point (geslon, geslat)
  
  
  ! INPUT
  ! geslon    longitude coordinate of the central point
  ! geslat    latitude  coordinate of the central point
  ! rad       the radius around the central point ?? is different from trkrinfo%rad ?
  
  
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
  !           Note: except if very high latitude, feature removed
  !           < 0 : In case of a regional grid, the distance between 
  !                   the edge and the central point is inferior to rad 
  !                 In that case, the edge of the subgrid is fixed to
  !                   the edge of the main one


    implicit none
  
    real,    intent(in)  :: geslon, geslat, rad
    integer, intent(inout) :: ibeg, jbeg, iend, jend, ret
    integer fix_i, fix_j, ip, jp, i, j
    real dist, degrees

    ret = 0
    
  !-----------------
  !-- to the East --
  !-----------------
    fix_i = int((geslon - glonmin + 2*dx)/dx)
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
          ! ret = 72
          exit
        else ! wrapping
          i = ip - imax
        endif
      endif
      
      call calc_dist_dir(geslon, geslat, glon(i), geslat, .false., &
                          dist, degrees)
                               
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
    fix_i = int((geslon - glonmin + dx)/dx)
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
          !ret = 72
          exit
        else ! wrapping
          i = ip + imax
        endif
      endif
      
      call calc_dist_dir(geslon, geslat, glon(i), geslat, .false., &
                          dist, degrees)
      
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
    fix_j = int((-geslat - glatmin + 2*dy)/dy)
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
      
      call calc_dist_dir(geslon, geslat, geslon, glat(j), .false., &
                          dist, degrees)
      
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
    fix_j = int((-geslat - glatmin + dy)/dy)
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
         
      call calc_dist_dir(geslon, geslat, geslon, glat(j), .false., &
                          dist, degrees)
      
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

  subroutine check_contour (ptx_i, ptx_j, fxy, cminmax, contour_min, &
                forced, mask_object, mask_object_cc, &
                contour_max, rcontour_max,ret)

! It tests if the contour "contour_min" around the point (ptx_i, ptx_j) 
!   is closed. The contour is not closed if it hits the edge of
!   a regional grid, if one of the point in it is already considered
!   as under the influence of another object or have a value under/above
!   the value at ptx_i, ptx_j.

! If closed, it returns the value of the largest closed contour
!   (contour_max), the area it includes (mask_object_cc), the area with
!   a constant gradient towards that latter (mask_object), and the mean
!   distance between the central point (ptx_i, ptx_j) and the largest
!   closed contour.

! Algorithm for cminmax = "min" ("max"): It first tests if the points
!   around ptx_i, ptx_j have a value below(above) contour_min and are
!   either hiden (by masked_out), have a value under the one of the
!   central point (ptx_i, ptx_j), or is outside the gird (condition for
!   no closed contour). If not, all the points are flagged as part of
!   the area under the influence of the object (mask_object). All their
!   neighbours which have not already been flagged are then considered 
!   for the no closed contour conditions. Then, if their value is 
!   below(above) contour_min or above(below) the value of at least one 
!   of their flagged neighbours, they are themselves flagged (flagging 
!   conditions). The algorithm continues to consider the next neighbours.
!   until one of the condition for no closed contour is met or until no
!   new points are being flagged (when therefore contour_min is a closed
!   contour)

! If forced = .true., the algorithm continues even if a condition for
!   no closed contour is found. At the end, contour_max is fixed to
!   contour_min


! Note: in that implementation, the central point needs to be on the
!    grid, it could be adapted for any points.
! Note: masked_out might be an argument if different than "area influenced by other object"
! previously check_closed_contour


  ! INPUT
  ! ptx_i     index coordinate along i of the point tested (longitude)
  ! ptx_j     index coordinate along j of the point tested (latitude)
  ! fxy       field considered (often data_detection)
  ! cminmax   either "max" or "min", indicates whether ptx is a min or a max
  ! forced    either .true. or .false., indicates whether to compute 
  !             mask_object, if a closed contour is not found
  ! contour_min   the value of the minimum contour to be tested
    
  ! OUPUT
  ! mask_object   if a contour_min is closed or forced = .true.,
  !                 includes (value = .true.) all the points with a  
  !                 value below(above) contour_min, not masked and  
  !                 connected to the central point. It also includes all
  !                 the points with a value above(below) contour_min
  !                 and  connected to one of the previous points by a
  !                 positive(negative) gradient.
  !               Especially usefull for first_ges_center since it  
  !                 doesn't includes any local minimum(maximum).
  ! contour_max   maximum(minimum) value above(below) contour_min
  !                 forming a closed contour around the central point 
  !                 and not including any local minimum(maximum)
  ! mask_object_cc    includes all the points within the contour
  !                     contour_max
  ! rcontour_max      mean distance between the central point and the 
  !                     edge of mask_object_cc
  ! ret   return success of the subroutine
  !         0:  contour_min is a closed contour (if forced = .true.)
  !         61: contour_min hits the edge of a regional grid
  !         62: contour_min includes a point hiden by another object
  !         63: contour_min includes a point with a value below/above
  !               the one at the center
  

    implicit none

    integer,     intent(in) :: ptx_i, ptx_j
    real,        intent(in) :: fxy(:,:), contour_min
    character*3, intent(in) :: cminmax
    logical,     intent(in) :: forced

    logical, intent(out) :: mask_object(:,:)
    logical, intent(out) :: mask_object_cc(:,:)
    real,    intent(out) :: contour_max, rcontour_max
    integer, intent(out) :: ret

    real :: xcentval
    real contour, dist, degrees
    integer ipt, pt_max, ptn, iptn, i, j, nb
    integer pt_i(imax*jmax), pt_j(imax*jmax), ptn_i(8), ptn_j(8)
    logical contour_close, on_border
    integer mm ! 1: "max", -1: "min"


    ret = 0
    xcentval = fxy(ptx_i,ptx_j)

    if ( verb .ge. 3 ) then
      print *,' '
      print *,'--------------------------------------------------------'
      print *,'* Top of check_contour, ix= ',ptx_i,' jx= ',ptx_j
      print *,'--------------------------------------------------------'
      print *,'cminmax: ', cminmax
      print *,'fxy(ix,jx)= ', xcentval
      print *,'minmial contour= ', contour_min
    endif
    



    ! to simplify the if statements
    if (cminmax == 'min') then
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

    ! This array will represent all the points bellow contour_min connected
    !   to the central point or with a negative gradient towards that later area
    mask_object = .false.
    mask_object(ptx_i, ptx_j) = .true.

    !nptb = 0 ! number of possible points at the border

    contour_close = .true.


  !-------------------------------------
  !-- Loop on the points (pt_i, pt_j) --
  !-------------------------------------
    scan_contour: do while (ipt <= pt_max)
      !if ( verb .ge. 3 ) then
      !  print *, " "
      !  print *, " the following point is beeing treated: "
      !  print *, " -lon: ", glon(pt_i(ipt))
      !  print *, " -lat: ", glat(pt_j(ipt))
      !endif 
    
    
      call get_neighbours(pt_i(ipt), pt_j(ipt), fileinfo%grid_type, &
              imax, jmax, ptn, ptn_i, ptn_j)
          
      ! we hit the border of a regional grid
      if (fileinfo%grid_type == 'regional' .and. ptn /= 8) then
        if ( (contour_min-fxy(pt_i(ipt), pt_j(ipt)))*mm < 0)  then
        
          if ( verb .ge. 3 ) then
            print *,'We reach the border of a regional grid'
            print *,'The minimum contour is not closed'
          endif
          
          contour_close = .false.
          
          if (.not. forced) then
            ret = 61
            exit scan_contour
          endif
          
        else
          ptn = 0 ! its neighbours are not considered
        endif
        
      endif          


    !-----------------------------------
    !-- scan the 8 neighboring points --
    !-----------------------------------
      scan_neighbour : do iptn=1,ptn
      

        if (masked_out(ptn_i(iptn), ptn_j(iptn))) then
          if ( (contour_min-fxy(pt_i(ipt), pt_j(ipt)))*mm < 0 ) then
          
            if ( verb .ge. 3 ) then
              print *,'We reach the border of another object'
              print *,'The minimum contour is not closed'
            endif
            
            contour_close = .false.

            if (.not. forced) then
              ret = 62
              exit scan_contour
            endif
            
          endif


        else if ((xcentval-fxy(ptn_i(iptn), ptn_j(iptn)))*mm < 0) then
          if ( (contour_min-fxy(pt_i(ipt), pt_j(ipt)))*mm < 0 ) then
          
            if ( verb .ge. 3 ) then
              print *, "We reach a point with a lower/higher value ", &
                        "than the central point"
              print *, "fxy(i,j) :", fxy(ptn_i(iptn), ptn_j(iptn))
              print *, "cminmax: ", cminmax
              print *,'The minimum contour is not closed'
            endif
            
            contour_close = .false.

            if (.not. forced) then
              ret = 63
              exit scan_contour
            endif

          endif
          
        ! There is nothing to do if the point is already in mask_object
        else if (.not. mask_object(ptn_i(iptn), ptn_j(iptn))) then
          if (( (contour_min-fxy(ptn_i(iptn), ptn_j(iptn)))*mm <= 0 ) .or. &
              ( (fxy(pt_i(ipt), pt_j(ipt)) - &
                 fxy(ptn_i(iptn), ptn_j(iptn)))*mm > 0 )) then
        
            mask_object(ptn_i(iptn), ptn_j(iptn)) = .true.
            pt_max = pt_max + 1
            pt_i(pt_max) = ptn_i(iptn)
            pt_j(pt_max) = ptn_j(iptn)
            
          endif
          
        endif
      enddo scan_neighbour

      ipt = ipt + 1
      
    enddo scan_contour



    if (.not. contour_close .and. .not. forced) then
    
      ! mask_object = .false. we keep this as none of the point already
      !   investigted could be the center of a new low: it will hit the 
      !   point investigated here, which has a lower value 
      !   (only in first_ges_center)
      mask_object_cc = .false.
      contour_max = xcentval
      rcontour_max = -999.
      return
    endif
    
  
  !------------------------------------------
  !-- determine the largest closed contour --
  !------------------------------------------
    if (contour_close) then 
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
      
    else ! forced and contour not closed
      contour_max = contour_min
    endif

    
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
          
          if (ptn == 8) then ! not at this edge of the grid
            do iptn=1,ptn
              if(.not. mask_object_cc(ptn_i(iptn),ptn_j(iptn))) then
                on_border = .true.
              endif
            enddo
          endif
          
          if(on_border) then
            call calc_dist_dir(glon(ptx_i), glat(ptx_j), glon(i), &
                    glat(j), .false., dist, degrees)
            nb = nb + 1
            rcontour_max = rcontour_max + dist
          endif
          
        endif
      enddo
    enddo

    rcontour_max = rcontour_max/nb


  end subroutine
  
  
  !---------------------------------------------------------------------------------
  !--------------- Subroutine searching for a minmax around a point-----------------
  !---------------------------------------------------------------------------------
  subroutine find_minmax(geslon, geslat, fxy, minmax, ptc_i, ptc_j, &
                clon, clat, cval, found)
  
! This subroutine looks for either local minimum or maximum of a field 
!     fxy around a point determine by its coordinates (geslon, geslat).
! The first step consists in reducing the area studied to a square
!     including all points within a certain distance of the point
!     (cf. get_ij_bounds).
! Then, all local minmax within that distance and not influenced by
!     another object are detected
! Finally, the minmax selected is the one maximising the formula 
!     (val-mean(A)/std(A))*((dist-rad)/rad)**2 where val is its value,
!     dist its distance to the central point, rad the maximum of this
!     distance and A the area of search.
! If no values are found, clon, clat and cval = -9999.
  
! Note: I didn't choose between the indices or the coordinate of the
!     center. I can add a bicubic interpolation which would give a 
!     better approximation

! previously find_maxmin


  ! INPUT
  ! geslon   longitude of the point around which the minmax is searched
  ! geslat   latitude  of the point around which the minmax is searched
  ! fxy      the field where to look for the minmax
  ! minmax   indicates whether to look for a minimum ("min") or a 
  !             maximum ("max")
  
  ! OUTPUT
  ! ptc_i   index along i of the minmax
  ! ptc_j   index along j of the minmax
  ! clon    longitude of the minmax
  ! clat    latitude  of the minmax
  ! cval    value of fxy at the minmax
  ! logical indicating if a point was found
  

    implicit none

    real, intent(in) :: geslon, geslat, fxy(:,:)
    character*3, intent(in) :: minmax
    
    real, intent(out) :: clon, clat
    integer, intent(out) :: ptc_i, ptc_j
    real, intent(out) :: cval
    logical, intent(out) :: found
    
    integer ibeg, jbeg, iend, jend ! the boundaries of the local grid
    real, allocatable :: fxy2(:,:) ! local version of fxy
    logical , allocatable :: mask_s(:,:) ! area of search
    integer iend2, imax2, jmax2, iw, ip ! for wrapping issue
    real dist,degrees, mean, sd
    integer nptm, ptm(2), iptm, ptm2(1) ! for the possible minmax
    integer, allocatable :: ptm_i(:), ptm_j(:) ! for the possible minmax
    integer ptn, ptn_i(8), ptn_j(8), iptn ! the neighbours
    real fact(imax*jmax)
    logical local_minmax
    integer i,j 

    integer err_a1, err_a2, err_a3, err_a4
    
    integer ret

    
    
    if ( verb .ge. 3 ) then
      print *, ""
      print *, "Looking for a minmax ..."
    endif
    
  !------------------------------------------------
  !-- Set up the mask where to look for a minmax --
  !-- It reduces the area studied -----------------
  !------------------------------------------------
    call get_ij_bounds(geslon, geslat, trkrinfo%rad, ibeg, jbeg, iend, jend, ret)
    
    if (ret == 71) then
      if ( verb .ge. 1 ) then
        print *, ""
        print *, "!!! ERROR 71 in find_minmax, from get_ij_bounds:"
        print *, "!!! The caracteristic size of the object provided by"
        print *, "!!! the user, rad (", trkrinfo%rad, ")"
        print *, "!!! is smaller than a gird point."
        print *, "!!! dx: ", dx, " dy: ", dy
        print *, "!!! No minmax can be found, tracking ends"
      endif
      STOP 71
    endif

    !if (ret == 72) then
    !  if ( verb .ge. 1 ) then
    !    print *, ""
    !    print *, "!!! ERROR 72 in find_minmax, from get_ij_bounds:"
    !    print *, "!!! The caracteristic size of the object provided by"
    !    print *, "!!! the user, rad (", trkrinfo%rad, ")"
    !    print *, "!!! is larger than 1/2 circonference of the Earth"
    !    print *, "!!! Too many minmax can be found, tracking ends"
    !  endif
    !  STOP 72
    !endif
          
    if ( verb .ge. 2 .and. ret < 0) then
      print *, ""
      print *, " WARNING in find_minmax, from get_ij_bounds:"
      print *, " The edge of a regional grid has been hit"
      print *, " This might lead to some problem in finding a minmax"
    endif
    

        
    if (ibeg > iend) then ! wrapping, Rq, ibeg > 0 by definition
    
      allocate(mask_s(iend-ibeg+1+imax, jend-jbeg+1),stat=err_a1)
      allocate(fxy2  (iend-ibeg+1+imax, jend-jbeg+1),stat=err_a2)
      if (err_a1 /= 0 .or. err_a2 /= 0) then 
        if (verb .ge. 1) then
          print *, ""
          print *, "!!! ERROR 73 in find_minmax: "
          print *, "!!! allocating mask_s or fxy2 failed"
          print *, "iend-ibeg+1+imax=", iend-ibeg+1+imax
          print *, "jend-jbeg+1=", jend-jbeg+1
        endif 
        STOP 73
      endif      
      
      fxy2((imax-ibeg+2):(iend-ibeg+1+imax),:) = fxy(1:iend, jbeg:jend)
      fxy2(1:(imax-ibeg+1),:) = fxy(ibeg:imax, jbeg:jend)
      
      iend2 = iend + imax ! for the next loop
      imax2 = iend-ibeg+1+imax ! for the call to get_neighbours
      jmax2 = jend-jbeg+1 ! for the call to get_neighbours
      
    else
    
      allocate(mask_s(iend-ibeg+1, jend-jbeg+1),stat=err_a1)
      allocate(fxy2(iend-ibeg+1, jend-jbeg+1),stat=err_a2)
      if (err_a1 /= 0 .or. err_a2 /= 0) then 
        if (verb .ge. 1) then
          print *, ""
          print *, "!!! ERROR 74 in find_minmax: "
          print *, "!!! allocating mask_s or fxy2 failed"
          print *, "iend-ibeg+1=", iend-ibeg+1
          print *, "jend-jbeg+1=", jend-jbeg+1
        endif 
        STOP 74
      endif  
            
      fxy2 = fxy(ibeg:iend, jbeg:jend)

      iend2 = iend ! for the next loop
      imax2 = iend-ibeg+1 ! for the call to get_neighbours
      jmax2 = jend-jbeg+1 ! for the call to get_neighbours
      
    endif
    
    
    ! fill mask_s, the mask of search, only trues will be analysed
    mask_s = .true.
    do i=ibeg,iend2 ! Rq: 0 < ibeg < iend2, wrapping or not
    
      !wrapping
      iw = i
      if (i > imax) iw = i - imax
      
      do j=jbeg,jend
      
        call calc_dist_dir(geslon, geslat, glon(iw), glat(j), .false., &
                            dist, degrees)
                            
        if (dist > trkrinfo%rad .or. masked_out(iw, j)) then
          mask_s(i - ibeg + 1, j - jbeg + 1) = .false.
        endif
      enddo
    enddo
    
    ! To compute the weights at the last step before mask_s changes
    mean = sum(fxy2, mask = mask_s)/count(mask_s)
    sd   = sqrt(sum((fxy2-mean)**2, mask = mask_s)/(count(mask_s)-1))
    
    !print *, ibeg, iend
    !print *, jbeg, jend
    


  !-----------------------------
  !-- search for local minmax --    
  !-----------------------------
    nptm = 0 ! number of minmax found
    allocate(ptm_i(imax2*jmax2),stat=err_a3)
    allocate(ptm_j(imax2*jmax2),stat=err_a4)
    if (err_a3 /= 0 .or. err_a4 /= 0) then 
      if (verb .ge. 1) then
        print *, ""
        print *, "!!! ERROR 75 in find_minmax: "
        print *, "!!! allocating ptm_i or ptm_j failed"
        print *, "imax2*jmax2=", imax2*jmax2
      endif 
      STOP 75
    endif      
        
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
      if (ptn < 8 .and. fileinfo%grid_type == 'regional') then
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
                fxy2(ptm(1), ptm(2))) then
              local_minmax = .false.
              exit
            endif
          else 
            if (fxy2(ptn_i(iptn), ptn_j(iptn)) > &
                fxy2(ptm(1), ptm(2))) then
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
    if (verb >= 3) then
      print *, "   lon    lat         val       dist   fact"
    endif
    fact = 0

    if (nptm > 0) then
      do iptm=1,nptm
      
        call calc_dist_dir(geslon, geslat, glon(ptm_i(iptm) + ibeg -1), &
                   glat(ptm_j(iptm) + jbeg - 1), .false., dist, degrees)
      
        fact(iptm) = ((fxy2(ptm_i(iptm), ptm_j(iptm))-mean)/sd)* &
                    ((dist-trkrinfo%rad)/trkrinfo%rad)**2
                    
                    
        if (verb >= 3) then
          print "(F7.2, F7.2, ES12.4, EN11.2, F7.3)", &
            glon(ptm_i(iptm) + ibeg -1), glat(ptm_j(iptm) + jbeg - 1), &
            fxy2(ptm_i(iptm), ptm_j(iptm)), dist, fact(iptm)
        endif
      enddo
      
      if (minmax == 'min') then
        ptm2 = minloc(fact(1:nptm))
      else 
        ptm2 = maxloc(fact(1:nptm))
      endif      
      
      ptc_i = ptm_i(ptm2(1)) + ibeg - 1
      ptc_j = ptm_j(ptm2(1)) + jbeg - 1
      
      clon = glon(ptc_i)
      clat = glat(ptc_j)
      cval = fxy2(ptm_i(ptm2(1)), ptm_j(ptm2(1)))

      if (ptc_i > imax) then  ! wrapping
        ptc_i = ptc_i - imax
      endif
      
      found = .true.
      
    else
      ptc_i = -999
      ptc_j = -999
    
      clon = -9999.
      clat = -9999.
      cval = -9999.
      
      found = .false.
    endif
    
    
     if ( verb .ge. 3 ) then
      if (found) then
        print *, "A new minmax has been found around (", &
                  geslon, geslat, ") at :"
        print *, "-lon: ", clon
        print *, "-lat: ", clat
        print *, "-val: ", cval        
      else
        print *, "No minmax found around (", geslon, geslat, ")"
      endif
    endif

    print *, "toto2"

  
  end subroutine
  


  
  !---------------------------------------------------------------------------------
  !----------------- Subroutine precising the center of an object ------------------
  !---------------------------------------------------------------------------------

  subroutine fixcenter (pt_i, pt_j, mask_object_max, &
                mask_object_min, mask_object, fix_lon, fix_lat, &
                fixvals, contour_max, rcontour_max)


! It precises the center of an object, depending on the minmax from the 
! fields data_fixcenter and recompute contour_max and rcontour_max if
! needed

! Note: This piece of code is especially useful to have a better estimate
!    where to find the next point of the track

  ! INPUT
  ! pt_i    index along i of the center found in data_detection
  ! pt_j    index along j of the center found in data_detection
  ! mask_object_max   output of find_minmax on data_detection (mask_object)
  ! mask_object_min   output of find_minmax on data_detection (mask_object_cc)
  
  ! OUTPUT
  ! fix_lon   longitude of the center
  ! fix_lat   latitude  of the center
  ! fixvals   values of the minmax found for the fields data_fixcenter
  ! mask_object   the actual mask of the object, to add to masked_out
  
  ! IN-OUT
  ! contour_max     equals contour_max from find_minmax on data_detection,
  !                   or minimum contour including all the centers from 
  !                   data_fixcenter (not closed in that case)
  ! rcontour_max   mean distance between the new center and the contour_max
  
  
    implicit none

    integer, intent(in) :: pt_i, pt_j ! the minmax from field detection
    logical, intent(in) :: mask_object_max(:,:), mask_object_min(:,:)
    
    real, intent(out) :: fix_lon,fix_lat,fixvals(:)
    logical, intent(out) :: mask_object(:,:)
    
    real, intent(inout) :: contour_max, rcontour_max
    
    
    real geslon, geslat ! the lon/lat of the center from the field detection
    
      ! the centers in the different fixcenter fields
    real clon(numfield%fixcenter), clat(numfield%fixcenter)
    integer ptc_i(numfield%fixcenter), ptc_j(numfield%fixcenter)
    logical found(numfield%fixcenter)
      
    real fix_lon1, fix_lat1 ! the center after simple average
    integer ifc, nfc ! number of fixcenter fields with a minmax
    integer nb, i, j, iptn, ptn, ptn_i(8), ptn_j(8), ccret
    real wt, wts, dist, degrees, val
    logical in_mask, on_border, mask_object_notused(imax, jmax)



    if ( verb >= 3 ) then
      print *, " "
      print *, "--------------------------------"
      print *, "* At the beginning of fixcenter "
      print *, "--------------------------------"
    endif
    
    
    geslon = glon(pt_i)
    geslat = glat(pt_j)
    
    
    if (.not. any(Rflag%fixcenter) .or. numfield%fixcenter == 0) then
      if ( verb >= 3 ) then
        print *, "No data are available for center fixing"    
        print *, "The center remains the same as previously detected:"
        print *, " -lon: ", geslon
        print *, " -lat: ", geslat
      endif
      fix_lon = geslon
      fix_lat = geslat
      mask_object = mask_object_min
      fixvals = -9999.
      return
    endif
      

     
  !-------------------------------------------------------------
  !-- First, we compute the minmax of the fields in fixcenter --
  !-------------------------------------------------------------
    do ifc = 1, numfield%fixcenter
      if (Rflag%fixcenter(ifc)) then
        call find_minmax(geslon, geslat, data_fixcenter(:,:,ifc), &
                fname%fc_minmax(ifc), ptc_i(ifc), ptc_j(ifc), &
                clon(ifc),clat(ifc),fixvals(ifc), found(ifc))
      else
        found(ifc) = .false.
      endif
    enddo
    
    
  !--------------------------------------------------------------
  !-- Print out a table listing of the locations of the minmax --
  !--------------------------------------------------------------
    if ( verb .ge. 3 ) then
      print *, ' '
      print *, '--------------------------------------------------'
      print *, " Fixes have been made for the object which"
      print *, " initial guess point is"
      print *, " -lon: ", geslon
      print *, " -lat: ", geslat
      print *, " The centers are: "

      do ifc = 1, numfield%fixcenter
        if(found(ifc)) then
            print *, " For ", trim(fname%fixcenter(ifc)), "-", &
                      fname%fixcenter_lev(ifc)
          print *, " lon: ",glon(ptc_i(ifc))," lat: ",glat(ptc_j(ifc)), &
                   " val: ", fixvals(ifc)
        endif
      enddo

    endif


  !--------------------------------
  !-- Calculate the new center --
  !--------------------------------
    ! First, a mean of all the centers are done
    ! Then, it's recalculated by giving more weight to centers closer 
    !   to the previously calculated mean position
    ! If some centers appeared to be at more than rad of the first mean
    !   position, they are not considered for the second one
    ! Note: we take a simple lon/lat average, without considering
    !   the Earth curvature
    ! A center will always be computed, even if it's the one provided


    ! simple average
    fix_lon = geslon
    fix_lat = geslat
    
    nfc = 0
    do ifc = 1, numfield%fixcenter
      if (found(ifc)) then
        nfc = nfc + 1
        fix_lat = fix_lat + glat(ptc_j(ifc))
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
    call calc_dist_dir(fix_lon1, fix_lat1, geslon, geslat, .false., &
            dist, degrees)
              
    if (dist > trkrinfo%rad) then ! priority to the detection field
      fix_lon1 = geslon
      fix_lat1 = geslat
      wt = 1
    else
      wt = ((dist-trkrinfo%rad)/trkrinfo%rad)**2
    endif
      
    wts = wt
    fix_lon = geslon * wt
    fix_lat = geslat * wt
      
    
    nfc = 0
    do ifc = 1, numfield%fixcenter
      if (found(ifc)) then
      
        call calc_dist_dir(fix_lon1, fix_lat1, glon(ptc_i(ifc)), &
                glat(ptc_j(ifc)), .false.,  dist, degrees)
        
        if (dist > trkrinfo%rad) then ! the parameter is not considered
          if (verb >= 3) then
            print *, " The variable ", trim(fname%fixcenter(ifc)), "-", &
                    fname%fixcenter_lev(ifc)
            print *, " appears to be too far from the others: ", dist
            print *, " it is discarded"
          endif
          fixvals(ifc) = -9999.
          clon(ifc) = -9999.
          clat(ifc) = -9999.
          
          exit
        endif
        
        wt = ((dist-trkrinfo%rad)/trkrinfo%rad)**2
        wts = wts + wt
        
        fix_lat = fix_lat + glat(ptc_j(ifc))*wt
        
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
    
    contour_max = data_detection(pt_i, pt_j)
    in_mask = .true.
    
    do ifc = 1,numfield%fixcenter
      if (found(ifc)) then
      
        if (.not. mask_object_min(ptc_i(ifc), ptc_j(ifc))) then
          val = data_detection(ptc_i(ifc), ptc_j(ifc))
        
          if (fname%detec_minmax == 'min') then
            contour_max = max(contour_max, val)
          else
            contour_max = min(contour_max, val)
          endif
        
          if(.not. mask_object_max(ptc_i(ifc), ptc_j(ifc))) then
            in_mask = .false.
          endif
          
        endif  
            
      endif
    enddo
    
    if(.not. in_mask) then
    
      call check_contour(pt_i, pt_j, data_detection, fname%detec_minmax, &
              contour_max, .true., mask_object_notused, mask_object, &
              contour_max, rcontour_max,ccret)
              
    else if (contour_max /= data_detection(pt_i, pt_j)) then
    
      mask_object = mask_object_max
      do i=1,imax
      do j=1,jmax
        if(mask_object(i,j) .and. data_detection(i,j) > contour_max) then
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
            
            if(on_border) then
              call calc_dist_dir(fix_lon, fix_lat, glon(i), glat(j), &
                                  .false., dist, degrees)
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
  !----- Subroutine that compute the movement speed and direction of an object -----
  !---------------------------------------------------------------------------------  
  subroutine get_speed_dir(iobj, its, fixlon, fixlat, mask_object, &
                              first_detection)  
  ! It fills speed and direction vectors defined in the module
  ! Note: speed is in m/s, as should be the advection field
  
    implicit none
    
    integer, intent(in) :: iobj, its
    real, intent(in) :: fixlon, fixlat
    logical, intent(in) :: first_detection, mask_object(:,:)
    
    real dist, degree ! in calcdist
    real extrapoled_direction, extrapoled_speed, dt
    logical extrapol ! true if the extrapolation was computed
    real wind_speed, wind_direction
    real wind_speed_ij, wind_u, wind_u_ij, wind_v, wind_v_ij
    integer n_ij ! number of selected points to compute the mean wind
    real wts, wt ! weightening for the mean wind
    integer i, j, iadv
    
    
    
    if (verb >= 3) then
      print *, ""
      print *, "Computing speed and direction ..."
    endif 
    
   
  !------------------------------------
  !-- Method 1: Linear extrapolation --
  !------------------------------------
  ! First, just do a simple linear extrapolation from the previous fix
  !   position through the current fix position. 
  ! Doesn't work for newly detected low

     
    if (.not. first_detection) then ! the object previous position is known
      if (prevlat(iobj) >= -90 .and. prevlon(iobj) >= 0) then ! To be sure they are initialised
      
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
        
        call calc_dist_dir(prevlon(iobj), prevlat(iobj), fixlon, &
                fixlat, .true., dist, degree)
                
        extrapoled_direction = degree
        extrapoled_speed = dist * 1000 / dt ! this is m/s
        extrapol = .true.
        
      else
        if (verb >= 1) then 
          print *, ""
          print *, "WARNING in get_speed_dir"
          print *, "the previous latitude or longitude of an object is"
          print *, "not initialised while it is not its first detection"
          print *, "prevlat: ", prevlat(iobj)
          print *, "prevlon: ", prevlon(iobj)
        endif
        extrapol = .false.
        
      endif
    else
      extrapol = .false.
    endif


    
  !------------------------------
  !-- Method 2: Wind advection --
  !------------------------------
  ! We take an average of the wind in the area influenced by the object
  !   and at less than rad from the center.
  ! More weights are given to the wind closer to the center.
  ! If more than one field is provided to compute the advection, we take
  !   an average of it
  
    wind_speed = 0
    wind_u = 0
    wind_v = 0
    n_ij = 0
    wts = 0
    
  
    do i=1,imax
      do j=1,jmax
        if (mask_object(i,j)) then
        
          call calc_dist_dir(glon(i), glat(j), fixlon, fixlat, .true., &
                dist, degree)
                
          if (dist < trkrinfo%rad) then

            wind_speed_ij = 0
            wind_u_ij = 0
            wind_v_ij = 0
            n_ij = n_ij + 1
            
            do iadv=1,numfield%advection
              if (Rflag%advection(iadv)) then
                wind_speed_ij = wind_speed_ij + &
                  sqrt(data_u(i,j,iadv)**2 + data_v(i,j,iadv)**2)
                wind_u_ij = wind_u_ij + data_u(i,j,iadv)
                wind_v_ij = wind_v_ij + data_v(i,j,iadv)
              endif
            enddo

            wt = ((dist-trkrinfo%rad)/trkrinfo%rad)**2
            wts = wts + wt
            wind_speed = wind_speed + &
                         wind_speed_ij/count(Rflag%advection)*wt
            wind_u = wind_u + wind_u_ij/count(Rflag%advection)*wt
            wind_v = wind_v + wind_v_ij/count(Rflag%advection)*wt
            

          endif
        endif
        
      enddo
    enddo
    
    wind_speed = wind_speed / wts
    ! no need tonormalised wind_u / wind_v

    ! 360 is the North, 0 means no wind
    if (wind_u == 0.0) then
      if (wind_v == 0.0) then 
        wind_direction = 0.0
      else if (wind_v > 0) then
        wind_direction = 360.0
      else if (wind_v < 0) then
        wind_direction = 180.0
      endif
    else if (wind_u > 0.0) then
      if (wind_v == 0.0) then 
        wind_direction = 90.0
      else 
        wind_direction = 90. - atan(wind_v/wind_u)  / dtr
      endif
    else if (wind_u < 0.0) then
      if (wind_v == 0.0) then 
        wind_direction = 270.0
      else
        wind_direction = 270. - atan(wind_v/wind_u) / dtr
      endif
    endif         
    
    
  !-----------------------------
  !-- average the two methods --
  !-----------------------------
  
    if(extrapol) then
      speed(iobj) = (wind_speed + extrapoled_speed)/2
      direction(iobj) = (wind_direction + extrapoled_direction)/2
      if (abs(wind_direction - extrapoled_direction) > 180) then ! across 0°
        direction(iobj) = direction(iobj) - 180
        if (direction(iobj) < 0) then
          direction(iobj) = direction(iobj) + 360
        endif
      endif
     else
      speed(iobj) = wind_speed
      direction(iobj) = wind_direction
    endif
    
    if ( verb .ge. 3 ) then
      print *, ""
      print *, "Speed and direction for the object: ", &
                names_obj(iobj)
      print *, "wind direction: ", wind_direction
      print *, "wind_speed: ", wind_speed
      if (extrapol) then
        print *, "extrapoled direction: ", extrapoled_direction
        print *, "extrapoled speed: ", extrapoled_speed
      else
        print *, "No extrapolation was computed"
      endif
      print *, "therefore,"
      print *, "object direction (in °): ", direction(iobj)
      print *, "object speed (in km/", trim(fileinfo%time_unit), &
                  "): ", speed(iobj)
    endif
                  
      

  end subroutine


  
  !---------------------------------------------------------------------------------
  !----- Subroutine searching for new local min / max in the "detection" field -----
  !---------------------------------------------------------------------------------
  subroutine first_ges_center (its, ilast_object, masked)

! This subroutine searches for new object within the boundary of the box
!   of search fixed by the user in tha namelist. The areas masked by
!   previous object (already under tracking) are not considered. If a
!   new objest if found, call the function to name it (newobject_name),
!   to fix its center (fixcenter) and print the outputs in file.
  
  
  ! INPUT
  ! its    the timestep considered
  ! masked collect from main subroutines all the mask_object_max
  !        computed in check_contour from the previously tracked object.
  !        There can not be any new minmax there
  
  ! INOUT
  ! ilast_object   the total number of object caracterised during the
  !                   run of the program
  
    use test
    implicit none
    
    integer, intent(in)    :: its
    logical, intent(in)    :: masked(:,:)
    integer, intent(inout) :: ilast_object
    
    integer ibeg, iend, jbeg, jend
    integer iptn, ptn, ptn_i(8), ptn_j(8) ! from get_neighbours
    logical mask_detec(imax, jmax), mask_object(imax,jmax)
    logical mask_object_max(imax,jmax), mask_object_min(imax,jmax)
    real contour, minmax_val, contour_max, rcontour_max
    real fixlon, fixlat, speed, direction
    real vals_fc(numfield%fixcenter), vals_int(numfield%intensity)
    integer ccret, ptx(2), ilast_object_old, iint
    integer ptc_i(numfield%intensity), ptc_j(numfield%intensity) ! not used
    real clon(numfield%intensity), clat(numfield%intensity) ! not used
    logical found(numfield%intensity) !not used    
    
    
    if ( verb .ge. 2 ) then
      print *,' '
      print *,'*-------------------------*'
      print *,'* Top of first_ges_center *'
      print *,'*-------------------------*'

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
      print *, "The indices of the boundary box are :"    
      print *, "jbeg: ", jbeg
      print *, "jend: ", jend
      print *, "ibeg: ", ibeg
      print *, "iend: ", iend
    endif


  !-----------------------------------------------------
  !-- define mask_detec, where to look for new minmax --
  !-----------------------------------------------------
    ! mask_detec hide the areas where the objects already being tracked are
    mask_detec = .not. masked

    ! hide the areas beyond the box of search
    if (ibeg < 1) then
      mask_detec((iend+1):(ibeg-1),:) = .false.
    else if (iend > imax) then
      if (ibeg > imax + 1) then ! This case is only possible when grid is global and glonmin<0 
        mask_detec(1:(ibeg-1-imax),:) = .false.
        mask_detec((iend+1-imax):imax,:) = .false. ! we suppose glonmax>0 and thus iend<2*imax
      else
        mask_detec((iend+1-imax):(ibeg-1),:) = .false.   
      endif
    else ! no wrapping case
      if (ibeg > 1) then
        mask_detec(1:(ibeg-1),:)    = .false.
      endif
      if (iend < imax) then
        mask_detec((iend+1):imax,:) = .false.
      endif
    endif
    
    if (jbeg > 1) then
      mask_detec(:,1:(jbeg-1))    = .false.
    endif
    if (jend < jmax) then
      mask_detec(:,(jend+1):jmax) = .false.
    endif

    ! hide the border of the domain
    mask_detec(:,1)    = .false.
    mask_detec(:,jmax) = .false.
    if (fileinfo%grid_type == 'regional') then
      mask_detec(1,:)    = .false.
      mask_detec(imax,:) = .false.
    endif

  ilast_object_old = ilast_object
  


  !---------------
  !-- Main loop --
  !---------------
    ! It exits only when the minmax value not masked are above/below the
    ! threshold, or when all the points are masked
    search_loop: do while (any(mask_detec))
  
      if (fname%detec_minmax == 'min') then
        ptx = minloc(data_detection, mask = mask_detec)
        minmax_val = data_detection(ptx(1), ptx(2))
     
        if (minmax_val > trkrinfo%detec_thresh) then
          exit search_loop
        endif
                
      else
        ptx = maxloc(data_detection, mask = mask_detec)
        minmax_val = data_detection(ptx(1), ptx(2))
        
        if (minmax_val < trkrinfo%detec_thresh) then
          exit search_loop
        endif

      endif

      if ( verb .ge. 3 ) then
        print *,""
        print *,""
        print *,'Found a possible max/min at ptx(1)= ', &
                  ptx(1),' ptx(2)= ',ptx(2)
      endif
      
    if (fname%detec_minmax == 'min') then
      contour = minmax_val + trkrinfo%contint
    else
      contour = minmax_val - trkrinfo%contint
    endif
          
      ! The point is considered if a closed contour is found around it
      call check_contour(ptx(1), ptx(2), data_detection, &
              fname%detec_minmax, contour, .false., mask_object_max, &
              mask_object_min, contour_max, rcontour_max,ccret)
              
              
      ! if no closed contour
      if (ccret /= 0) then
        if ( verb .ge. 3 ) then
          print *, "No contour has been found: no new object"
        endif        

        ! We don't want to find this local minimum and the points
        !   already detected as being under his influence
        mask_detec = mask_detec .and. .not. mask_object_max
          
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
          print *, "that timestep; ts_tim(its): ", print_ts_tim(its)
          print *, "Using the record of the last objects being"
          print *, "tracked will discard all the old storms"
        endif  
        STOP 51
      endif        


      ilast_object = ilast_object + 1
      if ( verb .ge. 3 ) then
        print *, "A contour has been found: this is a new object"
        print *, "contour_max: ", contour_max
        print *, "rcontour_max: ", rcontour_max
      endif


      ! better approximation of the center
      call fixcenter (ptx(1), ptx(2), mask_object_max, &
                mask_object_min, mask_object, fixlon, fixlat, vals_fc, &
                contour_max, rcontour_max)


      ! name object & save the center's coordinates
      call newobject_name(its, ilast_object, fixlon, fixlat)
      prevlon(ilast_object) = fixlon
      prevlat(ilast_object) = fixlat
      minmax_detec(ilast_object) = minmax_val
      beingTracked(ilast_object) = .true.
      


      ! compute the minmax for the fields intensity
      if (verb >= 3) then
        print *, ""
        print *, "Intensity minmax:"
        print *, "-----------------"
      endif
      do iint = 1, numfield%intensity
        if (Rflag%intensity(iint)) then
          call find_minmax(fixlon, fixlat, data_intensity(:,:,iint), &
                  fname%int_minmax(iint), ptc_i(iint), ptc_j(iint), &
                  clon(iint),clat(iint),vals_int(iint), found(iint))
        else
          found(iint) = .false.
          ptc_i(iint) = -999
          ptc_j(iint) = -999
          clon(iint)  = -9999.
          clat(iint)  = -9999.
          vals_int(iint) = -9999.
        endif
      enddo
      

      

      ! compute speed and direction from wind speed
      call get_speed_dir(ilast_object, its, fixlon, fixlat, &
                           mask_object, .true.)
         
      ! write ouputs
      call default_output(ilast_object, its, contour_max, &
                            rcontour_max, vals_fc, vals_int)


      ! fill mask_detec with new false and masked_out with true
      mask_detec = mask_detec .and. .not. mask_object_max
      masked_out = masked_out .or. mask_object


    enddo search_loop


    if ( verb .ge. 3 ) then
      print *, " "
      print *, " For the timestep ", its
      print *, " number of new objects discovered: ", ilast_object - &
                                                       ilast_object_old
    endif
  
  
  end subroutine
  

end module




