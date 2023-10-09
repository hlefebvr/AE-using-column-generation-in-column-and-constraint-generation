#!/bin/bash

if [ "$(whoami)" = "utr_lefebvre" ]
then
  PROJECT_DIRECTORY=/home/utr_lefebvre/min-max-min-nested-cg/code
else
  PROJECT_DIRECTORY=/home/henri/CLionProjects/min-max-min-nested-cg/code
fi

INSTANCE_DIRECTORY=data
BUILD_DIRECTORY=cmake-build-debug
EXECUTABLE=FLP/FLP_solve
EXPERIMENTS_DIRECTORY=experiments

COUNTER=0

for FILE in $PROJECT_DIRECTORY/$INSTANCE_DIRECTORY/*
do
  for GAMMA in 4
  do

    for STD_PHASE_TIME_LIMIT in 10800 60 120
    do

      for WITH_HEURISTIC in false true
      do

        ARGS="$PROJECT_DIRECTORY/$EXPERIMENTS_DIRECTORY/run_one.sh $PROJECT_DIRECTORY/$BUILD_DIRECTORY/$EXECUTABLE $FILE $STD_PHASE_TIME_LIMIT $GAMMA $WITH_HEURISTIC"

        echo "Submitting $ARGS"

        if [ "$(whoami)" = "utr_lefebvre" ]
        then
          sbatch $ARGS
        else
          $ARGS
        fi

        COUNTER=$(($COUNTER+1))

        if [ $STD_PHASE_TIME_LIMIT -eq 10800 ]
        then
          break
        fi

      done

    done
  done
done

echo "Submitted ${COUNTER}"
