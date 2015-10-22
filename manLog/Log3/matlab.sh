#!/bin/sh
#PBS -q serial
#PBS -l nodes=1:ppn=2
#PBS -l walltime=12:00:00
#PBS -N BAFBAtestcase

outputDir="/home/zhenyal/rProject/Output"
model="E coli iJR904"
task_fun="BAFBA"

cd $PBS_O_WORKDIR
module load matlab
matlab -nodesktop -nodisplay -nosplash -r "math_task('$outputDir','$model','$task_fun');exit" 
