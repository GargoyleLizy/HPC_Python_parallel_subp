import sys
import argparse
import subprocess
import time

from mpi4py import MPI
#from enum import Enum
from collections import namedtuple
import copy
#from collections import deque

from taskmaster import taskmaster

# a helper function to show the result
# the mode variance is used for unifying different methods
# the workers in "master slave" start from 1 
# while the wokrers in "bulk" start from 0.
def show_result_int(worker_num,task_list):
    workers=[]
    for idx in range(size):
        worker_rank = idx
        temp_worker = {'worker_rank':worker_rank,'total_work_time':float(0),'total_tasks':0}
        workers.append(dict( temp_worker ) )
    
    #print(workers)
    undone_tasks=[]
    for idx in range(len(task_list)):
        if(task_list[idx]['status']!=0):
            undone_tasks.append ( task_list[idx] )
        else:
            #print('time',task_master.task_list[idx]['time'])
            for widx in range(len(workers)):
                if(workers[widx]['worker_rank']==task_list[idx]['worker']):
                    workers[widx]['total_work_time']+= task_list[idx]['time']
                    workers[widx]['total_tasks']+=1

    print("Undone tasks: %d."%len(undone_tasks))
    for idx in range(len(undone_tasks)):
        print(undone_tasks[idx])

    print("%d workers status:"%len(workers)) 
    for widx in range(len(workers)):
        print("Worker %d did %d tasks with %.4f time."% ( workers[widx]['worker_rank']
            ,workers[widx]['total_tasks'],workers[widx]['total_work_time'] ) )

# A helper function to decide wether the external process spend too much time on a task.
# Source:
#   http://stackoverflow.com/questions/4158502/python-kill-or-terminate-subprocess-when-timeout
def subprocess_execute(command, time_out):
    """executing the command with a watchdog"""

    # launching the command
    c = subprocess.Popen(command)

    # now waiting for the command to complete
    t = 0
    while t < time_out and c.poll() is None:
        time.sleep(1)  # (comment 1)
        t += 1

    # there are two possibilities for the while to have stopped:
    if c.poll() is None:
        # in the case the process did not complete, we kill it
        c.terminate()
        # and fill the return code with some error value
        returncode = -1  # (comment 2)

    else:                 
        # in the case the process completed normally
        returncode = c.poll()

    return returncode 





# ***** Important Preparation here *****

# ------ Constant Variables ------
# ExtnFunction define which external function was called to process the task
# Might want to put that part into input. to simplify test work.
ExtnFunction="matlab"
#ExtnFunction="echo"

# Mode decide which mode is used for the tasks distribution and processing
# It can be specified by the input. 
# Currently, only 'ms' master-slave and 'bulk' distribute tasks as bulks 
# are legimate.
Mode='ms'
#Mode='bulk'

# Wall time is used for master-slave architecture
# After the last job is distributed out. The master would start a timer
# If some node did not complete their tasks within the wall time.
# Master would decide that they are crashed and mark the tasks as undone
# by default, Wall time is None, and the master would not check it.
Walltime = None


# Master Interval decide the interval time of master sleeps
# It is seconds
Master_Interval = 2

# tag used for mpi communication
#class Tags(Enum):
#    ready=1
#    start=2
#    done=3
#    exit=4
Tags=namedtuple('Tags',['ready','start','done','exit'])
tags=Tags(1,2,3,4)

#Task=namedtuple('Task',['idx','task'])
#Result=namedtuple('Result',['idx','return_value','execution_time'])


# ----- Initialize MPI variables -----
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
name = MPI.Get_processor_name()
status = MPI.Status()   


# ***** Main Function Here Starts *****

# ------ Handle Input -----
# Input arguments should contain:
# 1. -o An output directory to save results
# 2. -i An input file that contains the tasks
# 3. -p A prject Directory that specify where to find essential resources.
# 4. -ef The external function that called by subprocess
# 5. -m The tasks distributing mode.
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--ifile')
    parser.add_argument('-o','--odir')
    parser.add_argument('-p','--pdir')
    parser.add_argument('-ef','--exfun')
    parser.add_argument('-m','--mode')
    parser.add_argument('-wt','--walltime')
    args=parser.parse_args()
   
# If essential arguments are missing, warn user and exit
if(args.odir==None or args.ifile == None or args.pdir ==None):
    print("Essential output_dir/input_file/project_dir is missing")
    print("Program exit.")
    exit()

# If 'exfun' is specified, use that exfun instead of predifined matlab.
# The exfun could be "echo" for test
if(args.exfun != None):
    ExtnFunction=args.exfun

# Two structure 'ms' and 'bulk' are allowed
if(args.mode !=None):
    Mode=args.mode

# wall time set the upper time limit of a single task done by worker 
# It needs users being familiar with the tasks and can estimate upper limit for 
# processing a single task. 
# Be care with this var. If the cluster is reliable, just do not set it.
if(args.walltime != None):
    Walltime=float( args.walltime )


if(Mode == "ms"):
    # original source for mpi4py-master-slave architecture can be found in 
    #   https://github.com/jbornschein/mpi4py-examples
    # I modified the structure to fit in this project.
    
    # Recording time
    start_time = time.time()

    if(rank==0):
        # task_master is an abstract manager to handle the tasks.
        task_master = taskmaster(args.odir,args.ifile,args.pdir)
        task_master.extract_tasks()

        num_workers = size -1
        closed_workers = 0
        last_timer_time = None
        print("Master starting with %d workers" % num_workers) 
    
        while closed_workers < num_workers:
            # decide whether break loop jump to result without waiting unfinished tasks
            if comm.Iprobe(source=MPI.ANY_SOURCE,tag=MPI.ANY_TAG):
                data = comm.recv(source=MPI.ANY_SOURCE,tag=MPI.ANY_TAG,status=status)
                source = status.Get_source()
                tag = status.Get_tag()
                if(tag==tags.ready):
                    if(task_master.remain_task()):
                        temp_task = task_master.pop_task()
                        comm.send(temp_task,dest=source,tag=tags.start) 
                        print("[Master][Send] task [%d] to worker [%d] "% (temp_task[0],source ))
                        task_master.task_list[temp_task[0]]['worker']=source
                    else:
                        comm.send(None,dest=source,tag=tags.exit)
                elif(tag == tags.done):
                    result=data
                    print("[Master][Recv] task [%d] from worker [%d]" % (result[0],source))
                    task_master.task_list[result[0]]['status']=result[1]
                    task_master.task_list[result[0]]['time']=result[2]
                    task_master.task_list[result[0]]['worker']=source
                elif(tag==tags.exit):
                    print("[Master][Recv] worker [%d] exit" % source)
                    closed_workers +=1

            else:
                if Walltime != None:
                    if (not task_master.remain_task()):
                        if(last_timer_time != None):
                            if (time.time()-last_timer_time) > Walltime:
                                print("last timer time: ",time.time()-last_timer_time)
                                break;
                        else:
                            last_timer_time = time.time()
                #time.sleep(Master_Interval)
                #print("no sended packt")
        total_time = time.time() - start_time
        print("Master finishing, total_time %f. detailed results showing below:" %total_time)
        show_result_int(size-1,task_master.task_list)

    else:
        # worker processes 
        while(True):
            comm.send(None,dest=0,tag=tags.ready)
            task = comm.recv(source=0,tag=MPI.ANY_TAG,status=status)
            tag=status.Get_tag()
            if(tag==tags.start):
                
                node_start_time = time.time()
                print(task[1]) 
                if(Walltime !=None):
                    temp_result = subprocess_execute([ExtnFunction,task[1]],Walltime)
                else:
                    temp_result = subprocess.check_call([ExtnFunction,task[1] ])
                
                temp_execution_time = time.time()-node_start_time

                task_result = (task[0],temp_result,temp_execution_time)

                # when testing node crash, could uncomment the code below.
                # This is not exactly what crash would go, but I did not come up a better idea
                comm.send(task_result,dest=0,tag=tags.done)
            elif(tag==tags.exit):
                break
        comm.send(None,dest=0,tag=tags.exit)

# bulk mode is much simpler, each node does its own work and 
elif(Mode=="bulk" or Mode=="bulk-r"):
    if(rank==0):
        start_time = time.time()

    task_master = taskmaster(args.odir,args.ifile,args.pdir)
    task_master.extract_tasks()
    
    # shuffle the tasks
    if Mode == "bulk-r":
        task_master.shuffle_tasks()

    task_master.chunk_tasks(size)
    # each node executing the tasks here
    sub_task_list = task_master.get_idx_chunk(rank,ExtnFunction)
    sub_result_list = []
    for idx in range(len(sub_task_list)):
        temp_start_time = time.time()
        temp_task = sub_task_list[idx][2]
        temp_task_idx = sub_task_list[idx][3]
        if(Walltime==None):
            temp_result = subprocess.check_call([ExtnFunction,temp_task])
        else:
            temp_result = subprocess_execute([ExtnFunction,temp_task],Walltime)

        temp_execution_time = time.time()-temp_start_time
        task_master.task_list[temp_task_idx]['time'] = temp_execution_time
        task_master.task_list[temp_task_idx]['status'] = temp_result
        sub_result_list.append(task_master.task_list[temp_task_idx])



    
    # gather result
    result_list = comm.gather(sub_result_list,root=0)
    if(rank==0):
        total_time = time.time() - start_time
        final_list = []
        for idx in range(len(result_list)):
            final_list+=result_list[idx]
        
        print("Total time in rank 0 is %f sec(s)."%total_time)
        show_result_int(size,final_list)

else:
    print("Input Mode \"%s\" is wrong, only \'ms\' and \'bulk\' allowed."%Mode)
    

