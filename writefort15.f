      program writefort15

      USE netcdf

      implicit none

c     date format : hour since timeref, as provided in fort.11
c     timestep in hours too
      integer startdate,enddate,timestep
      integer ncid, dimid, length, ntime, i, time
      integer, allocatable :: values(:)
      integer ierror, ret, ret2, ret3, ret4





c     opening and reading the namelist
      open(unit = 15, file = 'namelist.15', status = 'old', 
     &      action = 'read', iostat = ierror)

      if (ierror /= 0) then

        print *,''
        print *,'Warning : error when oppening namelist.15'
        print *,'We will use all timesteps of the netcdf file'

        startdate = -1
        enddate   = -1
        timestep  = -1

      else

        read(15,*) startdate
        read(15,*) enddate
        read(15,*) timestep

      endif

      close(15)


c     opening and reading the netcdf data file
      ierror = nf90_open('fort.11',nf90_nowrite,ncid)
      if (ierror /= nf90_noerr) then
        print *,''
        print *,'Warning : Error when opening fort.11 with netcdf tools'

        if (startdate == -1 .or. enddate == -1 .or. timestep == -1) then
          print *,'!!! ERROR No default value can be computed'
          STOP 21
        endif
      endif

      ret = nf90_inq_dimid(ncid, 'time', dimid)
      if (ret == nf90_noerr) then
        ret2 = nf90_inquire_dimension(ncid, dimid, len = length)
        if (ret2 == nf90_noerr) then
          allocate(values(length), stat = ret3)
          if (ret3 == 0) then
            ret4 = nf90_get_var(ncid, dimid, values)
            if (ret4 == nf90_noerr) then
              go to 101
            endif
          endif
        endif
      endif

      print *,''
      print *,'!!! ERROR when reading fort.11'
      STOP 22

 101  continue

c     if no information in the namelist
      if (startdate == -1) then
        startdate = values(1)
      endif

      if (enddate == -1) then
        enddate = values(length)
      endif

      if (timestep == -1) then
        timestep = values(2) - values(1)
      endif


 
c     filling fort.15 file
      ntime = (enddate - startdate) / timestep + 1

      open(unit = 9, file = 'fort.15', status = 'new', 
     &     action = 'write', iostat = ierror)

      if (ierror /= 0) then
        print *,''
        print *,'!!! ERROR when oppening fort.15'
        STOP 23
      endif


      do i = 1,ntime

        time = (startdate + (i-1)*timestep)
        if (any(values == time)) then

          if (time*60 < 99999) then
            write(9,85) i, time*60
 85         format (i4,1x,i5)
          else
            print *,''
            print *,'!!! ERROR when filling fort.15'
            print *,'!!! one of the time value exceed the format'
            STOP 24
          endif

        else
          print *,''
          print *,'!!! ERROR when filling fort.15, '
          print *,'!!! one of the time value is not in the netcdf file'
          STOP 25
        endif


      enddo

      close(9)


      end





