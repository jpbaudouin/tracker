########################################################
########################################################
############ Post-Processing after tracking ############
########################################################
########################################################

# Created 06-10-2017

# It contains all the test for the function and class definition of WDTracking_package.R


###################################
############## TESTS ##############
###################################

track <- new("Track", id  = "test",
                      lat = (seq(0,20)/10)**2,
                      lon = seq(0,20) + 5,
                      intensity1 = sin(seq(0,20)/2),
                      intensity2 = abs(sin(seq(0,20)/2)),
                      dateStart  = as.POSIXct("2017-10-07 15:31:00" ,format = "%Y-%m-%d %H:%M:%S"),
                      timeStep   = 60*60*6)


track2 <- new("Track", id  = "test2",
                       lat = (seq(0,20)/10)**1.5,
                       lon = seq(0,20)/2 + 10,
                       intensity1 = (sin(seq(0,20)/2) + 1)**2,
                       intensity2 = (abs(sin(seq(0,20)/2)) + 1)**2,
                       dateStart  = as.POSIXct("2017-10-07 15:31:00" ,format = "%Y-%m-%d %H:%M:%S"),
                       timeStep   = 60*60*6)

validObject(track)

trackName(track)

trackFullPath(track)
trackPath(track)


# Not for user
trackDateStart(track)
trackTimeStep(track)
trackMaxInt1(track)
trackMaxInt2(track)
trackLength(track)
trackLines(track)


tracks <- new("ListTracks", tracks = list(track,track2))

tracksName(tracks)
tracksDateStarts(tracks)
tracksTimeStep(tracks)
tracksMaxInt1(tracks)
tracksMaxInt2(tracks)
tracksLength(tracks)
tracksPath(tracks)
tracksSpatialLines(tracks)

tracksSize(tracks)

toto<-readFort66("ERAI_tracking_1979-2007.txt")

limit <- Line(cbind(c(65, 65), c(20, 45)))
toto2 <- tracksSelect(toto, year = 2000, month = c(1,2) crossing = limit, propagation = "eastwards")
