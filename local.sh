#!/bin/bash
#---------------------------------------------------------
##########################################################
########## intialisation function for local run ##########
##########################################################
#---------------------------------------------------------


################
## PARAMETERS ##
################

member="ERAI" # "memb1 memb2 ..."
sdate=20000101 #19790101
chunksize=3 #348	# ~number of files
#chunknum=12
#launch_time=24 # (hours)
#pplaunch_time=1 # (hours)
#pplaunch_RAM=4 # (GB)
pplaunch_proc=1


####################
## Initialisation ##
####################

repository=$PWD
w_path=$(cat tracker.conf | grep WORKING_PATH | cut -d"=" -f2 )

if [ $w_path = "Default" ]
then
  echo "working_path is set to default :"
  name_exp=$(cat tracker.conf | grep NAME_EXP | cut -d" " -f3 )
  cd ../
  w_path=$PWD
  echo $w_path
fi

mkdir -p $w_path/$name_exp/exec

cd $w_path/$name_exp/exec


#######################
## Editing launch.sh ##
#######################

echo "#################"
echo "Editing launch.sh"
echo "#################"
echo "..."

cd $w_path/$name_exp/exec

sed -i "s#/esnas/autosubmit/a0ae/data#$repository#g" $repository/tracker.conf
#changing datapath for first run

for memb in $member
do
    i=1
    month=$(echo $sdate| cut -c5-6)

    start_year=$(( $(( $(echo $sdate| cut -c1-4) + $(( $(($i-1))*$chunksize )) /12 )) ))
    start_month=$(printf %02d $(($(( ${month#0} + $(( $(($i-1))*$chunksize )) %12 )) )) )
    start_day=$(echo $sdate| cut -c7-8)
    start_hour=00
    start_date=$start_year$start_month$start_day


    cp $repository/launch.sh launch_${memb}.sh

    # Replace the value from the file tracker.conf
    for var in $(cat launch_${memb}.sh | grep %*% | cut -d"%" -f2,4 )
    do


      line=$(cat $repository/tracker.conf | grep $var )
  
      if [ ${#line} -gt 0 ]
      then
        val=$(echo $line | cut -d" " -f3 )
        sed -i "s#%$var%#$val#g" launch_${memb}.sh
      fi

    done

    # Replace the value from this script
    sed -i "s#repository=.*#repository=$repository#g" launch_${memb}.sh
    sed -i "s#working_path=Default#working_path=$w_path#g" launch_${memb}.sh
    sed -i "s/%MEMBER%/${memb}/g" launch_${memb}.sh
    sed -i "s/%Chunk_START_DATE%/$start_date/g" launch_${memb}.sh
    sed -i "s/%Chunk_START_YEAR%/$start_year/g" launch_${memb}.sh
    sed -i "s/%Chunk_START_MONTH%/$start_month/g" launch_${memb}.sh
    sed -i "s/%Chunk_START_DAY%/$start_day/g" launch_${memb}.sh
    sed -i "s/%Chunk_START_HOUR%/$start_hour/g" launch_${memb}.sh
    sed -i "s/%CHUNKSIZE%/$chunksize/g" launch_${memb}.sh
    sed -i "s/module/#module/g" launch_${memb}.sh


    chmod 755 ./launch_${memb}.sh

done


#########################
## Editing pplaunch.sh ##
#########################

echo "###################"
echo "Editing pplaunch.sh"
echo "###################"
echo "..."

for memb in $member
do
  
  cp $repository/pplaunch.sh pplaunch_${memb}.sh
  # Replace the value from the file tracker.conf
  for var in $(cat pplaunch_${memb}.sh | grep %*% | cut -d"%" -f2,4 )
  do

    line=$(cat $repository/tracker.conf | grep $var )

    if [ ${#line} -gt 0 ]
    then
      val=$(echo $line | cut -d" " -f3 )
      sed -i "s#%$var%#$val#g" pplaunch_${memb}.sh
    fi



  done

  # Replace the value from this script
  sed -i "s#working_path=Default#working_path=$w_path#g" pplaunch_${memb}.sh
  sed -i "s#repository=.*#repository=$repository/R#g" pplaunch_${memb}.sh
  sed -i "s/%MEMBER%/$member/g" pplaunch_${memb}.sh
  sed -i "s/%NUMPROC%/$pplaunch_proc/g" pplaunch_${memb}.sh
  sed -i "s/module/#module/g" pplaunch_${memb}.sh


  chmod 755 pplaunch_${memb}.sh

done


##########################
## Parallel Jobs launch ##
##########################

echo "####################"
echo "Parallel Jobs launch"
echo "####################"
echo "..."

for memb in $member
do

  ./launch_${memb}.sh
  #./pplaunch_${memb}.sh

done
  





