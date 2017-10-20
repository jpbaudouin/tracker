########################################################
########################################################
############ Post-Processing after tracking ############
########################################################
########################################################

# Created 06-10-2017

# It contains all the function and class definition for WDTracking.R


#################################################################################################
#################################################################################################

#require(sp)
require(stringr)
require(rgeos) # remove gIntersect for an sp function?

## Class definition
#####################

setClass("Track",
         slots = c(id  = "character",
                   lat = "numeric",
                   lon = "numeric",
                   intensity1 = "numeric",
                   intensity2 = "numeric",
                   dateStart  = "POSIXct",
                   timeStep   = "numeric"), # in second
         validity = function(object) {
                      isValid <- TRUE
                      if (min(length(object@lat), length(object@lon),
                              length(object@intensity1), length(object@intensity2)) !=
                          max(length(object@lat), length(object@lon),
                              length(object@intensity1), length(object@intensity2))) {
                        isValid <- FALSE
                        print("Error on the length of the lat, lon or ints slots")
                      }
                      if (max(object@lat) > 90 | min(object@lat) < -90 | min(object@lon) < 0) {
                        isValid <- FALSE
                        print("Error on the value of lat or lon")
                      }
                      return(isValid)
                    })


setClass("ListTracks",
         slots = c(tracks  = "list"),
         validity = function(object) {
                      isValid <- TRUE
                      names <- c()
                      # Test class is Track
                      for (track in object@tracks) {
                        if (class(track) != "Track") {
                        isValid = FALSE
                        break
                        }
                        names <- c(names,trackName(track))
                      }
                      # Test all names are unique
                      if (length(unique(names)) != length(object@tracks)) {
                        isValid = FALSE
                      }
                      return(isValid)
                    })

## Accessors for Track
########################

trackName <- function(object) object@id

trackFullPath <- function(object) {
                   data.frame(lat = object@lat,
                              lon = object@lon,
                              intensity1 = object@intensity1,
                              intensity2 = object@intensity2,
                              date = object@dateStart + ((1:length(object@lat)-1)*object@timeStep))}

trackPath <- function(object) {
               data.frame(lat = object@lat,
                          lon = object@lon)}


# Not for user
trackDateStart <- function(object) object@dateStart
trackTimeStep  <- function(object) object@timeStep
trackMaxInt1   <- function(object) max(object@intensity1)
trackMaxInt2   <- function(object) max(object@intensity2)
trackLength    <- function(object) length(object@intensity2)
trackLines     <- function(object) Lines(Line(cbind(object@lon, object@lat)), ID=object@id)


## Integrated Accessors for ListTracks
########################################

tracksName       <- function(object) sapply(object@tracks, trackName)
tracksDateStarts <- function(object) lapply(object@tracks, trackDateStart)
tracksTimeStep   <- function(object) sapply(object@tracks, trackTimeStep)
tracksMaxInt1    <- function(object) sapply(object@tracks, trackMaxInt1)
tracksMaxInt2    <- function(object) sapply(object@tracks, trackMaxInt2)
tracksLength     <- function(object) sapply(object@tracks, trackLength) #Length of each track
tracksPath       <- function(object) lapply(object@tracks, trackPath)
tracksSpatialLines <- function(object)  SpatialLines(lapply(object@tracks, trackLines))

tracksFullPath   <- function(object) {
                      list <- lapply(object@tracks, trackFullPath)
                      names(list) <- tracksName(object)
                      return(list)
                    }

tracksSize <- function(object) length(object@tracks)



## Replacement methods
########################
#`trackName<-` <- function(object, value) {
#  object@id = value
#  object
#}

tracksSel <- function(object, sel) {
  if (max(sel) > length(object@tracks)) print("error in the selection, out of indices")
  object@tracks <- object@tracks[sel]
  return(new("ListTracks", tracks = object@tracks))
}


## Read GFDL tracker fort.66 format file
##########################################
readFort66 <- function(file) {

  rawResult <- read.table(file, na.strings = c("   -999", " -999", "  -9999"), sep = ",")

  hour1 <- as.numeric(as.character(rawResult$V7))
  testh <- sort(unique(hour1))
  time_step <- min(testh[-1] - testh[-length(testh)])*3600 # Need test time steps don't change
  
  day <- as.Date(trunc(hour1 / 24), origin = as.Date(substr(rawResult$V4, 1, 8), format = "%Y%m%d"))
  hour <- hour1 %% 24
  
  date = as.POSIXct(paste(day," ", hour, ":00:00", sep = ""),format = "%Y-%m-%d %H:%M:%S")
  
  lat <- as.numeric(str_sub(rawResult$V8, -4, -2)) / 10
  lat[which(str_sub(rawResult$V8,-1,-1) == "S")] <- -lat[which(str_sub(rawResult$V8,-1,-1) == "S")]
  
  lon <- as.numeric(str_sub(rawResult$V9, -5, -2)) / 10
  lon[which(str_sub(rawResult$V9,-1,-1) == "W")] <- -lon[which(str_sub(rawResult$V9,-1,-1) == "W")] + 360


  name <- gsub(" ","",as.character(rawResult$V3))
  unique_name <- unique(name)
  #result <- data.frame( name, date, lat, lon, Z = rawResult$V11, vo = rawResult$V28)

  list_track = c()
  for (i in 1:length(unique_name)) {

    if (i%%100 == 0) {
      cat("\rtrack", i, "/",length(unique_name))
      flush.console()
    }

    track <- which(name == unique_name[i])

    # To fix the crossing of the longitude boundary : all longitude are > 0 and continuous
    lon_track <- lon[track]
    lon_diff <- lon_track[2:length(lon_track)] - lon_track[2:length(lon_track)-1]
    jump <- which(abs(lon_diff) > 180)
    for (j in jump) {
      if (lon_diff[j] < 0) {
        lon_track[(j+1):length(lon_track)] <- lon_track[(j+1):length(lon_track)] + 360
      } else {
        lon_track[(j+1):length(lon_track)] <- lon_track[(j+1):length(lon_track)] - 360
      }
    }

    if (min(lon_track) < 0) {
      lon_track <- lon_track + 360
    }

    list_track <- append(list_track,
                         new("Track", id  = unique_name[i],
                                      lat = lat[track],
                                      lon = lon_track,
                                      intensity1 = rawResult$V11[track],
                                      intensity2 = rawResult$V28[track],
                                      dateStart  = date[track][1],
                                      timeStep   = time_step))
  }

  return(new("ListTracks", tracks = list_track))

}


## Selection of tracks 
########################
tracksSelect <- function(object, length = 1, year = NA, month = NA, crossing = NA, propagation = "none") {

# length is in hour
# crossing is a Line(cbind(lon, lat))
  

  # time length
  sel_time <- which(tracksTimeStep(object)*(tracksLength(object)-1)/3600 >= length)

  object <- tracksSel(object,sel_time)

  # year
  if (is.numeric(year)) {
  
    sel_year <- which(sapply(tracksDateStarts(object),format, format = "%Y") %in% year)

    object <- tracksSel(object,sel_year)
  }

  # month
  if (is.numeric(month[1])) {
    sel_month <- which(sapply(tracksDateStarts(object),format, format = "%m") %in% str_pad(month, 2, pad = "0"))

    object <- tracksSel(object,sel_month)
  }

  # Test eastwards propagation : lon_beg<lon_end
  if (propagation == "eastwards") {
    test_eastward <- function(dataFrame) {
      if(dataFrame$lon[1] < tail(dataFrame$lon,1)) return(TRUE)
      else return(FALSE)
    }

    sel_eastwards <- which(sapply(tracksPath(object), test_eastward) == TRUE)

    object <- tracksSel(object,sel_eastwards)

  }

  # Test Crossing the longitude 65 between 20 and 45
  if (class(crossing) == "Line") {
    crossing2 <- crossing
    crossing2@coords[,1] <- crossing2@coords[,1] + 360

    sel_crossing <- gIntersects(tracksSpatialLines(object),
                                 SpatialLines(list(Lines(crossing, ID="a"))),
                                 byid = TRUE) |
                    gIntersects(tracksSpatialLines(object), 
                                 SpatialLines(list(Lines(crossing2, ID="a"))),
                                 byid = TRUE)

   sel_crossing <- which(sel_crossing == TRUE)

    object <- tracksSel(object,sel_crossing)
  }

  return(object)

}

