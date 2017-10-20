#!/bin/bash
#----------------------------------------------------------
###########################################################################
########## First draft of the post-processing launching function ##########
###########################################################################
#----------------------------------------------------------



###########################
##### Loading modules #####
###########################
module purge
module load PROJ
module load GEOS
module load GDAL
module load R


######################
##### Parameters #####
######################
# Coming from tracker.conf

working_path=%WORKING_PATH%
member=%MEMBER%

# R parameters
name_exp=%NAME_EXP%
scale=%SCALE%
lat_max=%LAT_MAX%
time_min=%TIME_MIN%
wind_min=%WIND_MIN%
max_landfall=%MAX_LANDFALL%
no_cores=%NUMPROC%
user=%moore_USER%

# Others
repository=/esnas/autosubmit/%EXPID%/proj/BSC-Cyclone-Tracker-tools/R


if [ $working_path = "Default" ]
then
  echo "working_path is set to default :"
  working_path=/esnas/scratch/Earth/$user/Tracker
  echo $working_path
fi


mkdir -p $working_path/$name_exp/$member
cd $working_path/$name_exp/$member

# Making the results file
rm -f result.txt
for file in $(ls ../$member.*/result*.txt)
do
  cat $file >> result.txt
done

if [ ! -f result.txt ]
then
  echo "no file result.txt : exit"
  exit
fi

mkdir -p plot/hurricane
mkdir -p file

ln -s $repository/base_fun.R base_fun.R
  
$repository/resultProcessing.R $name_exp  $scale $lat_max $time_min $wind_min $max_landfall $no_cores

if [ ! -f plot/*season.eps ]
then
  echo "error of the launching of the post-processing"
  exit
fi

rm base_fun.R
rm result.txt
rm -r -f ../$member.*



