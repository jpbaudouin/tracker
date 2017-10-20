#!/usr/bin/env Rscript

########################################################
########################################################
############ Post-Processing after tracking ############
########################################################
########################################################

# Created 05-10-2017


source("WDTracking_package.R")

require(sp)
require(ggplot2)


#source("R/base_fun.R")
  
#require(stringr)
#require(maps)
#require(maptools)
#require(raster)
#require(rgeos)
#require(rworldmap)
#require(geosphere)
#require(grid)
#require(gridExtra)
#require(gtable)S
#require(abind)


######################
##### Parameters #####
######################





#################################################################################################
#################################################################################################
##### CORE SCRIPT 
###################


######################
## Reading the data ##
######################

ERAI_tracks <-readFort66("ERAI_tracking_1979-2007.txt")

ERAI_tracks <-readFort66("ERAI_tracking3_2000.txt")
  
########################## 
## Selecting the storms ##
##########################

limit <- Line(cbind(c(70, 70), c(20, 45)))
test <- tracksSelect(ERAI_tracks, length = 24, year = 2000, month = 1, crossing = limit, propagation = "none")

tracksFullPath(test)

###################################################################################################
## Plot
########

result5 <- array(NA, dim = c(4, length(names_select), max(t)))
for (n in 1:length(names_select)) {
  data <- result[which(result$name == names_select[n] &), ]
  result5[1, n, 1:dim(data)[1]] = data$lat
  result5[2, n, 1:dim(data)[1]] = data$lon
  result5[3, n, 1:dim(data)[1]] = data$Z
  result5[3, n, 1:dim(data)[1]] = data$vo
}

  
p <- plot_track(data = result5[c(1,2,4), , , drop = FALSE], type = "wind", limits = c(-180, 180, 0, 90), title = "Tracks detected", legend = "vorticity", return_plot = TRUE)










