#!/bin/bash

if [ "$(whoami)" = "utr_lefebvre" ]
then
  PROJECT_DIRECTORY=/home/utr_lefebvre/AE-using-column-generation-in-column-and-constraint-generation/code/
else
  PROJECT_DIRECTORY=/home/henri/CLionProjects/min-max-min-nested-cg/code
fi

INSTANCE_DIRECTORY=JSP/data
BUILD_DIRECTORY=cmake-build-debug
EXECUTABLE=JSP/JSP_solve
EXPERIMENTS_DIRECTORY=JSP/experiments

COUNTER=0

for FILE in $PROJECT_DIRECTORY/$INSTANCE_DIRECTORY/*
do
  for GAMMA in 5 7 9
  do

    for STD_PHASE_TIME_LIMIT in 0 10800 # 120 60
    do

      for WITH_HEURISTIC in false true
      do

        for NON_OPTIMAL_PRICING in false
        do

          ARGS="$PROJECT_DIRECTORY/$EXPERIMENTS_DIRECTORY/run_one.sh $PROJECT_DIRECTORY/$BUILD_DIRECTORY/$EXECUTABLE $FILE $STD_PHASE_TIME_LIMIT $GAMMA $WITH_HEURISTIC $NON_OPTIMAL_PRICING"

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

        if [ $STD_PHASE_TIME_LIMIT -eq 10800 ]
        then
          break
        fi

      done

    done
  done
done

echo "Submitted ${COUNTER}"
