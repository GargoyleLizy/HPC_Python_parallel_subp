import sys
import argparse
import subprocess
import time

#from collections import namedtuple
#from collections import deque


# ------ Task_Master ------ 
# A class responsible for fetching and distributing tasks 
# output_dir: All output of tasks would goto there.
# input_file_path: The input tasks should be put there
# 
class Task_Master:
    output_dir=""
    input_file_path=""
    project_path=""  
    task_list=[]
    task_current_indx=0
     
    # constant string pre
    matlab_argums_pre = "-nodesktop -nodisplay -nosplash -r "
    
    def __init__(self,output_dir,input_file_path,project_path):
        self.output_dir=output_dir
        self.input_file_path=input_file_path
        self.project_path=project_path

    def set_input(self,input_file):
        self.input_file_path=input_file
    
    #---Functions extract tasks from inputfile---
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
            elif(task_parts[1]=='BHFBA' or task_parts[1]=='DBFBA'):
                temp = {'model':task_parts[0], 'task_fun':task_parts[1],'rxnlist':task_parts[2]
                        ,'targetRxn':task_parts[3],'substrateRxn':task_parts[4]
                        ,'MaxKOs':int(task_parts[5]),'imax':int(task_parts[6])
                        ,'status':"undone",'time':0}
            self.task_list.append(temp)
        else:
            print "task is not correct formated:",task_str
    
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

    # --- Check remaining task number 
    def remain_task(self):
        return len( self.task_list ) > self.task_current_indx 


    
# ***** Main Function Here Starts *****

# Recording time
start_time = time.time()

# ------ Handle Input -----
# Input arguments should contain:
# 1. An output directory to save results
# 2. An input file that contains the tasks
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--ifile')
    parser.add_argument('-o','--odir')
    parser.add_argument('-p','--pdir')
    args=parser.parse_args()
    print "Input file : ", args.ifile
    print "Output dir : ",args.odir
    print "Project dir :",args.pdir
    task_master = Task_Master(args.odir,args.ifile,args.pdir)


# Debug: Just for Debug. clean these code later
with open(task_master.input_file_path,'r') as infile:
        for line in infile:
            print line.split('-')

# task_master extract tasks from input file
task_master.extract_tasks()

# Iterate through the task list and process them
while( task_master.remain_task() ):
    temp_start_time = time.time()
    
    matlab_arguments = task_master.pop_task()
    #print "!test: ", matlab_arguments
    if( subprocess.check_call(["matlab",matlab_arguments])  ):
        #print "failed?"
        task_master.task_list[task_master.task_current_indx-1]['status'] = "failed"
    else:
        #print "done?"
        task_master.task_list[task_master.task_current_indx-1]['status'] = "done"
    
    temp_execution_time = time.time() - temp_start_time
    task_master.task_list[task_master.task_current_indx-1]['time']=temp_execution_time

for task in task_master.task_list:
    print task

total_time = time.time()-start_time
print "total_time: ",total_time
