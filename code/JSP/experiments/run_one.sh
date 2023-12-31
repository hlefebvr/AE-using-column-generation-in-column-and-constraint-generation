#!/bin/bash
#SBATCH --nodes 1
#SBATCH --mem 80GB
#SBATCH --time 0-02:30:00
#SBATCH --mail-type FAIL
#SBATCH --ntasks 1 # 1 processor to be used
#SBATCH --constraint XEON_SP_6126

if [ "$(whoami)" = "utr_lefebvre" ]
then
  module purge
  module load gcc
fi

$@
