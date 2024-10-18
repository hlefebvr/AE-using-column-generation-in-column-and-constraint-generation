#!/bin/bash

if [ "$(whoami)" = "utr_lefebvre" ]
then
  PROJECT_DIRECTORY=/home/utr_lefebvre/AE-using-column-generation-in-column-and-constraint-generation/code/
else
  PROJECT_DIRECTORY=/home/henri/Research/AE-using-column-generation-in-column-and-constraint-generation/code/
fi

INSTANCE_DIRECTORY=JSP/data
BUILD_DIRECTORY=cmake-build-debug
EXECUTABLE=JSP/JSP_solve
EXPERIMENTS_DIRECTORY=JSP/experiments

COUNTER=0

for FILE in $PROJECT_DIRECTORY/$INSTANCE_DIRECTORY/*
do
  for GAMMA in 3 4
  do

    for STD_PHASE_TIME_LIMIT in 0 18000 # 120 60
    do

      for WITH_HEURISTIC in true # false true
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

          if [ $STD_PHASE_TIME_LIMIT -eq 18000 ]
          then
            break
          fi

        done

        if [ $STD_PHASE_TIME_LIMIT -eq 18000 ]
        then
          break
        fi

    done
  done
done

echo "Submitted ${COUNTER}"
