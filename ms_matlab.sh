#!/bin/sh
#PBS -q serial
#PBS -l nodes=1:ppn=4
#PBS -l walltime=16:00:00
#PBS -N ms01-N1C4-tasklist

outputDir="/home/zhenyal/rProject/Output/"
inputfile="/home/zhenyal/rProject/Input/tasklist"
projectDir="/home/zhenyal/rProject/"
#model="E coli iJR904"
#task_fun="BAFBA"
#rxn_list='all'
#targetRxn="UMPK"
#substrateRxn="UNK3"
#MaxKOs=5

cd $PBS_O_WORKDIR
module load matlab
module load python/3.2.3-gcc
module load openmpi-gcc
#matlab -nodesktop -nodisplay -nosplash -r "math_task('$outputDir','$model','$task_fun','$rxn_list','$targetRxn','$substrateRxn',$MaxKOs);exit" 
mpirun -np 4 python3 intP.py -o $outputDir -i $inputfile -p $projectDir -ef matlab -m ms
