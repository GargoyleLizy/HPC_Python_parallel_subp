from mpi4py import MPI
import re
import datetime
import sys

#------ initialize the mpi variables ------
comm = MPI.COMM_WORLD;
size = comm.Get_size();
rank = comm.Get_rank();
name = MPI.Get_processor_name();

# -- Recording time ---
if rank ==0:
    start_time = datetime.datetime.now()
    print("start_time: ",start_time)

# ---- handle the input ----
# First argument is input directory 
# Second  is output directory
# Third should be the external function that called
# Fourth should be the arguments provided to external function?
#   -> no, this should be handled by a seperate function.
if(len(sys.argv)>=3)
