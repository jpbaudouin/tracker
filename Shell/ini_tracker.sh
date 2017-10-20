#!/bin/bash
#--------------------------------
#################################
########## Compilation ##########
#################################
#--------------------------------


#########################
## Module & Paramaters ##
#########################

module load netCDF-Fortran/4.2-foss-2015a

repository=/esnas/autosubmit/%EXPID%/proj/BSC-Cyclone-Tracker-tools



###########################
## Compiling tracker.exe ##
###########################
echo " "
echo " "
echo "#######################"
echo "## Compiling library ##"
echo "#######################"
echo " "
echo " "

cd $repository/tracker_util
./clean -a  # from GFDL
cp configure-gnu.tracker_util  configure.tracker_util  # special configure file for gnu compiler not provided by GFDL
./compile   # from GFDL


echo " "
echo " "
echo "###########################"
echo "## Compiling tracker.exe ##"
echo "###########################"
echo " "
echo " "


cd $repository/gfdl-vortextracker

mkdir -p trk_exec

./clean -a  # from GFDL

cp configure-gnu.trk  configure.trk  # special configure file for gnu compiler not provided by GFDL
sed -i "s#LIB_W3_PATH    =   #LIB_W3_PATH    =   "$repository/tracker_util/libs"#g" configure.trk
sed -i "s#LIB_BACIO_PATH =   #LIB_BACIO_PATH =   "$repository/tracker_util/libs"#g" configure.trk
./compile   # from GFDL

if [ ! -f $repository/gfdl-vortextracker/trk_exec/tracker.exe ]
then
  echo "tracker.exe not compile"
  exit
fi


#####################
## writefort15.exe ##
#####################
echo " "
echo " "
echo "###############################"
echo "## Compiling writefort15.exe ##"
echo "###############################"
echo " "
echo " "

cd $repository

gfortran -I/shared/earth/software/netCDF-Fortran/4.2-foss-2015a/include writefort15.f -L/shared/earth/software/netCDF-Fortran/4.2-foss-2015a/lib -lnetcdff -o writefort15.exe


if [ ! -f $repository/writefort15.exe ]
then
  echo "writefort15.exe not compile"
  exit
fi

