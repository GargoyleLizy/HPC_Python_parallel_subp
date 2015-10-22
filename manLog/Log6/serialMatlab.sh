#!/bin/sh
#PBS -q serial
#PBS -l nodes=1:ppn=2
#PBS -l walltime=15:00:00
#PBS -N mathtest

outputDir="/home/zhenyal/rProject/Output/"
inputfile="/home/zhenyal/rProject/Input/tasklist1"
projectDir="/home/zhenyal/rProject/"
#model="E coli iJR904"
#task_fun="BAFBA"
#rxn_list='all'
#targetRxn="UMPK"
#substrateRxn="UNK3"
#MaxKOs=5

cd $PBS_O_WORKDIR
module load matlab
module load python/2.7.10-gcc
#matlab -nodesktop -nodisplay -nosplash -r "math_task('$outputDir','$model','$task_fun','$rxn_list','$targetRxn','$substrateRxn',$MaxKOs);exit" 
python serial.py -o $outputDir -i $inputfile -p $projectDir
