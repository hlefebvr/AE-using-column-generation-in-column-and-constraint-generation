#!/bin/bash
#SBATCH --nodes 1
#SBATCH --mem 32GB
#SBATCH --time 0-03:30:00
#SBATCH --mail-type NONE
#SBATCH --ntasks 1 # 1 processor to be used
#SBATCH --constraint XEON_SP_6126

if [ "$(whoami)" = "utr_lefebvre" ]
then
  module purge
  module load gcc
fi

echo "cmd: ${@}"

$@
