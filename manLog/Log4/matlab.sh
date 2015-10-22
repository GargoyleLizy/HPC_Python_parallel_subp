#!/bin/sh
#PBS -q serial
#PBS -l nodes=1:ppn=2
#PBS -l walltime=2:00:00
#PBS -N BAFBAtest

outputDir="/home/zhenyal/rProject/Output/"
model="E coli iJR904"
task_fun="BAFBA"
rxn_list='all'
targetRxn="UMPK"
substrateRxn="UNK3"
MaxKOs=5

cd $PBS_O_WORKDIR
module load matlab
module load python/2.7.10-gcc
#matlab -nodesktop -nodisplay -nosplash -r "math_task('$outputDir','$model','$task_fun','$rxn_list','$targetRxn','$substrateRxn',$MaxKOs);exit" 
python subp.py -o $outputDir
