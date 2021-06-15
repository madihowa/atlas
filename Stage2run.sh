#!/usr/bin/env bash

interactive=
DATE=$(date '+%F')
RUN=-1
counter=1
RUNPYTHON="1"
echo $date

while [ "$1" != "" ]; do
  case $1 in
    -d | --date )   shift
    DATE=$1
    ;;
    -r | --run )    shift
    RUN=$1
    ;;
  esac
  shift
done
echo "$date"

for FILE in `ls -l`
do
  if test -d $FILE
  then
    echo $FILE
    if [[ "$FILE" == *"$DATE"* ]]
    then
      if [[ "$FILE" == *run_"$RUN"* ]]
      then
          for FILE2 in `ls $FILE -l`
          do
            if [[ "$FILE2" == *root* ]]
            then
              RUNPYTHON=0
              echo "Found Root FILE not rerunning network"
              rm $FILE/$FILE2 #Remove root file then regenerate it later.
            fi
          done
      fi
      ((counter++))
    fi

  fi
done
#printf "\033c"
clear
# echo $RUNPYTHON

if (($RUNPYTHON == "1"));
then
  # echo $RUN
  python MasterStage2.py -r $RUN -d $DATE
fi

if (("$RUN" == "-1")); then
  #echo "New Run"
  RUN=$counter
fi

StringRes=Results_
StringRun=_run_


DIRECTORY="$StringRes$DATE$StringRun$RUN"
if (($RUNPYTHON == "1"));
then
  mv callback_history.csv $DIRECTORY
fi
#echo $DIRECTORY
#echo $DATE


cd $DIRECTORY
root -q ../runTStage2.C
cd ..
