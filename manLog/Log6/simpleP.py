import sys
import argparse
import subprocess
import time

from mpi4py import MPI
#from collections import namedtuple
#from collections import deque

# ------ Important Constant Variable ------
#ExtnFunction="matlab"
ExtnFunction="echo"

# ----- Initialize MPI variables -----
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
name = MPI.Get_processor_name()


# ------ Task_Master ------ 
# A class responsible for fetching and distributing tasks 
# output_dir: All output of tasks would goto there.
# input_file_path: The input tasks should be put there
# project_path : The path of essential resources for matlab function 
class Task_Master:
    output_dir=""
    input_file_path=""
    project_path=""  
    # stores tasks
    task_list=[]

    # Used to record tasks for better workload balancing
    task_current_indx=0

    # Used for simple workload distribution
    task_chunks=[]
    
    # constant string pre
    matlab_argums_pre = "-nodesktop -nodisplay -nosplash -r "
    
    def __init__(self,output_dir,input_file_path,project_path):
        self.output_dir=output_dir
        self.input_file_path=input_file_path
        self.project_path=project_path

    def set_input(self,input_file):
        self.input_file_path=input_file
    
    #---Functions extract tasks from inputfile---
    # This function assume each line represent a task
    # Modify it if this assumption does not fit your project.
    def extract_tasks(self):
        with open(self.input_file_path,'r') as infile:
            for line in infile:
                self.append_task(line)

    #--- !!! Modify append_task() to fit different projects
    def append_task(self,task_str):
        task_parts = task_str.split('-')
        if(len(task_parts)>=6):
            #temp = self.Task(task_parts[0],task_parts[1],task_parts[2],task_parts[3],task_parts[4],int(task_parts[5]),"undone",0)
            if(task_parts[1]=='BAFBA'):
                temp = {'model':task_parts[0], 'task_fun':task_parts[1],'rxnlist':task_parts[2]
                        ,'targetRxn':task_parts[3],'substrateRxn':task_parts[4]
                        ,'MaxKOs':int(task_parts[5]),'imax':0
                        ,'status':"undone",'time':0}
                self.task_list.append(temp)
            elif(task_parts[1]=='BHFBA' or task_parts[1]=='DBFBA'):
                temp = {'model':task_parts[0], 'task_fun':task_parts[1],'rxnlist':task_parts[2]
                        ,'targetRxn':task_parts[3],'substrateRxn':task_parts[4]
                        ,'MaxKOs':int(task_parts[5]),'imax':int(task_parts[6])
                        ,'status':"undone",'time':0}
                self.task_list.append(temp)
        else:
            print("task is not correct formated:",task_str)
    
    # --- !!! Modify pop_task() to fit different projects 
    # popout the rightmost task and cat them together to be called by subprocess.
    def pop_task(self):
        if(self.remain_task()):
            temp_task = self.task_list[self.task_current_indx]
            self.task_current_indx+=1
            matlab_argums_cont = "\"math_task("+ "\'"+self.project_path+"\'" \
                        +",\'"+ self.output_dir + "\'" \
                        +",\'"+ temp_task['model'] + "\'" \
                        +",\'"+ temp_task['task_fun'] +"\'" \
                        +",\'"+ temp_task['rxnlist'] + "\'" \
                        +",\'"+temp_task['targetRxn'] + "\'" \
                        +",\'" + temp_task['substrateRxn']+"\'" \
                        +"," + str(temp_task['MaxKOs']) \
                        +"," + str(temp_task['imax']) \
                        + ");exit\""

            return self.matlab_argums_pre + matlab_argums_cont
        else: 
            return "Done"
    
    # --- Distribute the tasks as blocks ---
    # This is a not so intelligent way to distribute the tasks
    def chunk_tasks(self,num_worker):
        self.task_chunks=[]
        avg = int( len(self.task_list)/num_worker)
        last=0
        for idx in range(num_worker):
            last+=avg
            self.task_chunks.append(last)
        self.task_chunks[num_worker-1]=len(self.task_list)
   
    def get_idx_task(self,task_idx):
        if(task_idx<len(self.task_list)):
            temp_task = self.task_list[task_idx]
            matlab_argums_cont = "\"math_task("+ "\'"+self.project_path+"\'" \
                        +",\'"+ self.output_dir + "\'" \
                        +",\'"+ temp_task['model'] + "\'" \
                        +",\'"+ temp_task['task_fun'] +"\'" \
                        +",\'"+ temp_task['rxnlist'] + "\'" \
                        +",\'"+temp_task['targetRxn'] + "\'" \
                        +",\'" + temp_task['substrateRxn']+"\'" \
                        +"," + str(temp_task['MaxKOs']) \
                        +"," + str(temp_task['imax']) \
                        + ");exit\""
            return self.matlab_argums_pre + matlab_argums_cont
        else: 
            return None

    
    # Execute the tasks with rank as input
    def execute_taskchunks(self,worker_idx):
        if(worker_idx==0):
            pre_idx=0
        else:
            pre_idx=self.task_chunks[worker_idx-1]
        
        for idx in range(pre_idx,self.task_chunks[worker_idx]):
            temp_start_time = time.time()
            
            subprocess.check_call([ExtnFunction,self.get_idx_task(idx)])
            
            temp_execution_time = time.time()-temp_start_time
            task_master.task_list[idx]['time']=temp_execution_time



    # --- Check remaining task number 
    def remain_task(self):
        return len( self.task_list ) > self.task_current_indx
    
    # --- Return total number of tasks
    def total_num_task(self):
        return len(self.task_list)

    
# ***** Main Function Here Starts *****

# Recording time
start_time = time.time()

# ------ Handle Input -----
# Input arguments should contain:
# 1. An output directory to save results
# 2. An input file that contains the tasks
# 3. A prject Directory that specify where to find essential resources.
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--ifile')
    parser.add_argument('-o','--odir')
    parser.add_argument('-p','--pdir')
    args=parser.parse_args()
    #print "Input file : ", args.ifile
    #print "Output dir : ",args.odir
    #print "Project dir :",args.pdir  



task_master = Task_Master(args.odir,args.ifile,args.pdir)
#task_master extract tasks from input file
task_master.extract_tasks()
task_master.chunk_tasks(size)
#print "size,",size, ";task list,", len(task_master.task_list)

#print "rank,",rank, ";task_chunk,",task_master.task_chunks[rank]
node_start_time = time.time()

task_master.execute_taskchunks(rank)

node_process_time = time.time() - node_start_time

print("[Total time] Node ",rank," takes ",node_process_time, " sec(s)")


# Iterate through the task list and process them
#while( task_master.remain_task() ):
#    temp_start_time = time.time()
#    
#    matlab_arguments = task_master.pop_task()
#    #print "!test: ", matlab_arguments
#    if( subprocess.check_call([ExtnFunction,matlab_arguments])  ):
#        #print "failed?"
#        task_master.task_list[task_master.task_current_indx-1]['status'] = "failed"
#    else:
#        #print "done?"
#        task_master.task_list[task_master.task_current_indx-1]['status'] = "done"
#    
#    temp_execution_time = time.time() - temp_start_time
#    task_master.task_list[task_master.task_current_indx-1]['time']=temp_execution_time

#for task in task_master.task_list:
    #print task


