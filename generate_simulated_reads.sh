#!/usr/bin/bash
#SBATCH --cpus-per-task 32 --time=12:00:00 --mem=64g
set -e

WD=/gpfs/gsfs12/users/lorenziha/KAREN_USDIN/TK_119/USDIN_LAB_TK_32/simulated_genome

time ${WD}/simuscop/bin/simuReads ${WD}/config_LB2BAV7596_wes.txt > simuReads.log 2>&1

echo =====================
echo Done!!
echo =====================
