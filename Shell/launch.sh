#!/bin/bash
#----------------------------------------------------------
###########################################################
########## First draft of the launching function ##########
###########################################################
#----------------------------------------------------------


###########################
##### Loading modules #####
###########################

module purge
module load netCDF-Fortran/4.2-foss-2015a
module load CDO
module load GEOS


################################
################################
########## Parameters ##########
################################
################################

# Coming from tracker.conf

type=%TYPE%
datapath=%DATAPATH%
working_path=%WORKING_PATH%
format=%FORMAT%
keepformateddata=%KEEP_FORMATTED_DATA%
name_exp=%NAME_EXP%

sdate15=%SDATE15%
edate15=%EDATE15%
tstep15=%TSTEP15%

wbd=%WESTBD%
ebd=%EASTBD%
nbd=%NORTHBD%
sbd=%SOUTHBD%

mslpth=%MSLPTHRESH%
mslpthgrad=%MSLPTHRESHGRAD%
v850=%V850THRESH%
contin=%CONTINT%
wc=%WCORE_DEPTH%
verb=%VERB%



# coming from autosubmit
repository=/esnas/autosubmit/%EXPID%/proj/BSC-Cyclone-Tracker-tools
member=%MEMBER%
sdate=%Chunk_START_DATE%
year=%Chunk_START_YEAR%
month=%Chunk_START_MONTH%
month_length=%CHUNKSIZE%
day=%Chunk_START_DAY%
hour=%Chunk_START_HOUR%
user=%moore_USER%


# Testing the parameters
###

echo "Testing the parameters..."

if [ $working_path = "Default" ]
then
  echo "working_path is set to default :"
  working_path=/esnas/scratch/Earth/$user/Tracker
  echo $working_path
fi

if [[ ! "netcdf grib" =~ $type ]]
then
  echo "Error on parameters type :" $type
  echo "should be either netcdf or grib"
  return 41
fi

if [[ ! "CMOR onefile rawgrib" =~ $format ]]
then
  echo "Error on parameters format :" $format
  echo "should be CMOR, onefile or rawgrib"
  return 43
fi

if [[ ! "T F" =~ $keepformateddata ]]
then
  echo "Error on parameters keepformateddata :" $keepformateddata
  echo "should be either T or F"
  return 44
fi

echo "Parameters OK"
echo " "



##############################
##############################
##### Launching function #####
##############################
##############################


launch () { # need an argument with the value :
            # 0 -> Cyclogenesis mode without input of previous cyclones
            # 1 -> Cyclogenesis mode WITH input of previous cyclones
            # 2 -> Tracking only mode (Thus with input of previous cyclones)

            # the second argument must be the name of the member to be analysed

  month=$(printf %02d $month)
  memb=$2
  echo " "
  echo "Treatment of time " $year$month
  echo "Launching mode" $1
  echo "Member " $memb

   

  ##############
  # FORMATTING #
  ##############

  ## CMOR
  if [ $format = "CMOR" ]
  then

    echo "Launching formatting of CMORised data"

    $repository/formatting/CMOR.sh $datapath $month $year $memb
    test=$?

    if [ $test -gt 1 ]
    then
      echo "CMORisation formatting has failed : " $test
      return 11
    fi

    mkdir -p ../data/
    mv data/$year$month.nc ../data/$memb"_"$year$month.nc
    file=../data/$memb"_"$year$month.nc

    echo "Formatting done"

  fi 

  ## GRIB
  if [ $format = "rawgrib" ]
  then

    echo "Launching formatting of raw grib"

    $repository/formatting/rawgrib.sh $datapath/$memb/outputs $month $year
    test=$?

    if [ $test -gt 1 ]
    then
      echo "formatting from raw grib has failed : " $test
      return 12
    fi
  
    mkdir -p ../data/
    mv data/$year$month.nc ../data/$memb"_"$year$month.nc
    file=../data/$memb"_"$year$month.nc

    echo "Formatting done"

  fi 

  ## NO FORMATTING
  if [ $format = "onefile" ]
  then

    if [ ! -f $datapath/$memb*$year$month.nc ]
    then
      echo "No file detected to proccess :"
      echo $datapath/$memb*$year$month.nc
      if [ $1 -eq 2 ]
      then
        echo "but it was for the tracking only mode (no error)"
        return
      else
        return 13
      fi
    fi
 
    file=$datapath/*$memb*$year$month.nc

    echo "Formatting not needed"

  fi 

  echo " "
  echo "####################### "
  echo "End of the formating"
  echo "File to be treated : " 
  echo $file
  echo "####################### "
  echo " "
    
  # Cleaning and making the files needed by the tracker
  #####################################################

  echo " "
  echo "cleaning"

  cd $working_path/$name_exp/$memb.$sdate/
  rm -f fort.14

  if [ $1 -eq 1 ]
  then

    if [ -f fort.67 ]
    then
      cp fort.67 fort.14
    else
      echo "No file fort.67 where must be the previously detected cyclones"
    fi

  elif [ $1 -eq 2 ]
  then

    if [ -f fort.67 ]
    then
      cp fort.67 fort.14
    else
      echo "No file fort.67 where must be the previously detected cyclones"
      echo "Running in tracking only mode not possible"
      return 30
    fi

  fi


  echo "initialisation"

  rm -f fort.11 fort.15 fort.66
  ln -s $file fort.11

  ./writefort15.exe

  if [ ! -f fort.15 ]
  then
    echo "error from writefort15.exe, no fort.15 file"
    return 20
  fi

  freq=$[$[$(sed -n 2p fort.15 | cut -c5-11) - $(sed -n 1p fort.15 | cut -c5-11)]/60*100]
  line=$(echo $(wc -l fort.15) | cut -c1-4)
  echo $line" timesteps for a frequency of " $freq" centahour have been found for fort.15"


  # making namelist
  cp $repository/namelist .
  sed -i "s/inp%bcc=/inp%bcc="$(echo $year| cut -c1-2)"/g" namelist
  sed -i "s/inp%byy=/inp%byy="$(echo $year| cut -c3-4)"/g" namelist
  sed -i "s/inp%bmm=/inp%bmm="$month"/g" namelist
  sed -i "s/inp%bdd=/inp%bdd=$day/g" namelist
  sed -i "s/inp%bhh=/inp%bhh=$hour/g" namelist
  sed -i "s/inp%filetype=/inp%filetype='$type'/g" namelist
  sed -i "s/atcfname=/atcfname='$name_exp'/g" namelist
  sed -i "s/atcfymdh=/atcfymdh="$year$month"0106/g" namelist
  sed -i "s/atcffreq=/atcffreq="$freq"00/g" namelist
  sed -i "s/trkrinfo%westbd=/trkrinfo%westbd=$wbd/g" namelist
  sed -i "s/trkrinfo%eastbd=/trkrinfo%eastbd=$ebd/g" namelist
  sed -i "s/trkrinfo%northbd=/trkrinfo%northbd=$nbd/g" namelist
  sed -i "s/trkrinfo%southbd=/trkrinfo%southbd=$sbd/g" namelist

  if [ $1 -eq 2 ]
  then
    sed -i "s/trkrinfo%type=/trkrinfo%type='tracker'/g" namelist
  else 
    sed -i "s/trkrinfo%type=/trkrinfo%type='tcgen'/g" namelist
  fi

  sed -i "s/trkrinfo%mslpthresh=/trkrinfo%mslpthresh=$mslpthgrad/g" namelist
  sed -i "s/trkrinfo%mslpthresh2=/trkrinfo%mslpthresh2=$mslpth/g" namelist
  sed -i "s/trkrinfo%v850thresh=/trkrinfo%v850thresh=$v850/g" namelist
  sed -i "s/trkrinfo%contint=/trkrinfo%contint=$contin/g" namelist
  sed -i "s/wcore_depth=/wcore_depth=$wc/g" namelist
  sed -i "s/verb=/verb=$verb/g" namelist

  # Launching tracker & aftermath
  ################################
  echo " "
  echo " "
  echo "LAUNCHING"
  echo " "
  echo " "

  ./tracker.exe < namelist

  if [ ! $? -eq 0 ]
  then
    echo "error after launching the tracker : exit"
    exit
  fi

  echo " "
  echo " "
  echo " "
  echo " "

  cat fort.66 >> result$sdate.txt

  if [ $keepformateddata = 'F' ] && [ $file = "../data/$year$month.nc" ]
  then
    rm ../data/$year$month.nc
  fi

  return 0

}



#########################################
#########################################
########## loop to call launch ##########
#########################################
#########################################


#########################################
##### Making the folder for the exp #####
#########################################

echo "Making the folder and the coping the file"

mkdir -p $working_path/$name_exp/$member.$sdate

cd $working_path/$name_exp/$member.$sdate/

rm -f writefort15.exe writefort15.f tracker.exe namelist namelist.15

ln -s $repository/writefort15.exe .
ln -s $repository/gfdl-vortextracker/trk_exec/tracker.exe .

echo $sdate15 > namelist.15
echo $edate15 >> namelist.15
echo $tstep15 >> namelist.15





## parameters

month=$[${month#0}-1]
rm -f result$sdate.txt


## Main loop over all the month of the jobs

for i in $(seq 1 $month_length)
do

  month=$[${month#0}+1]

  if [ $month -eq 13 ]
  then
    month=$[${month#0}-12]
    year=$[$year+1]
  fi

  launch 1 $member

  ret_launch=$?

  if [ $ret_launch -gt 1 ]
  then 
    echo "treatment of the data for " $year$month " has failed " $ret_launch
  fi

done


# Try to launch the following month in tracking-only mode

if [ $ret_launch -eq 0 ]
then
  month=$[${month#0}+1]

  if [ $month -eq 13 ]
  then
    month=1
    year=$[$year+1]
  fi

  launch 2 $member
fi








