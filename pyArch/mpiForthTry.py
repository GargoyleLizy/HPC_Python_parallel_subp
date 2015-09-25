from mpi4py import MPI
from collections import Counter
import sys
import re
import datetime

#---- initialize the mpi variables ---------------
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
name = MPI.Get_processor_name()

mode = MPI.MODE_RDONLY

#---- recording the time -------------------------
if rank == 0:
    start_time = datetime.datetime.now()
    # TIME RECODING
    print("start_time: ",start_time)

#---- initialize the file search variables -------
file_path = "data/Twitter.csv"

# if there is an command argument, take it as query word
# and the default query word is ash
query_word = "ash" 
if(len(sys.argv)>=2):
    if(sys.argv[1]):
        query_word = sys.argv[1]

# FOR DEBUGING
#print(rank,"(th)proc searching for:",query_word)

query_word_count = 0
user_counter = Counter()
topic_counter = Counter()

#---- sub functions ------------------------------

# chunk a list into certain num chunks with evenly size
def chunkIt(data_list, num_chunk):
    avg = len(data_list)/float(num_chunk)
    data_chunks = []
    last = 0.0

    while last < len(data_list):
        data_chunks.append(data_list[int(last):int(last+avg)])
        last += avg

    return data_chunks

# extract the twitter text from a csv line. by using the regular expression
# return a list - matchObj, and the text is matchObj(2).
def extractText( csvLine ):
   reFlags = re.M|re.I|re.U
   pattern = "\d+,\d+,\d+,\d+,(.*)text\"\":\"\"(.*)\"\",\"\"in_reply_to_status_id\"\":(.*)"
   matchObj = re.match(pattern, csvLine, reFlags)
   return matchObj

# count the frequence of query word in each line
# and return number.
def countTargetPerLine(query,twit):
    i = 0
    words = ''.join(c if c.isalnum() else ' ' for c in twit.lower()).split()

    for word in words:
        if(word == query):
            i = i+1
    return i

# find the @username in each line and return them as a list
def findUserPerLine(twit):
   users = re.findall(r"(?<=@)\w+",twit.lower())
   return users

# find the #topic in each line and return them as a list
def findTopicPerLine(twit):
   topics = re.findall(r"(?<=#)\w+",twit.lower())
   return topics

# add a list of terms to the Counter()
def addTermsToCounter(term_count,terms):
   line_count = Counter(terms)
   term_count = term_count + line_count
   return term_count

# count the query_word frequency, frequencies of usernames, frequencies
# of topics. return them as a list - [int,Counter(),Counter()]
def countAll(query, twit_list):
    query_count = 0
    user_count = Counter()
    topic_count = Counter()
    i = 0
   
    for twit in twit_list:
        i+=1
        if(i%100 == 0):
            print(i)
        twit_exrt = extractText(twit)
        if twit_exrt:
            twit_text = twit_exrt.group(2)
            
            query_count += countTargetPerLine(query,twit_text)

            users = findUserPerLine(twit_text)
            user_count = addTermsToCounter(user_count,users)
   
            topics = findTopicPerLine(twit_text)
            topic_count = addTermsToCounter(topic_count,topics)
   
    return [query_count, user_count, topic_count]

# get a list of lists. add each element of sublist together to get a single
# list
def listsToOne(result_lists,total_result):
    for result in result_lists:
        for i in range(0,len(result)):
            total_result[i] += result[i]
            
# print the result and adding favor text
def printResult(total_result,query_word):
    print("The word \"", query_word, "\" appeared ", total_result[0], "time(s) in the file")
    print("Top 10 users: ",total_result[1].most_common(10))
    print("Top 10 topics: ",total_result[2].most_common(10) )
    

#---- main mpi code here -------------------------
#---- rank0 get the file content -----------------

if rank == 0:
    twits = open(file_path,'r')
    data = twits.readlines()
    twits.close()
    data_chunks = chunkIt(data,size)
    IO_end_time = datetime.datetime.now()
    # TIME RECODING
    print("IO_end_time: ", IO_end_time)
else:
    data_chunks = None
    
# SCATTER data : each proc get its chunk of data
data_chunk = comm.scatter(data_chunks,root=0)

# For DEBUGING
#print(len(data_chunk))

proc_start_time = datetime.datetime.now()
print("proc %d start processing at %s" %(rank,proc_start_time))

# each node process its chunk
proc_result = countAll(query_word,data_chunk)

proc_end_time = datetime.datetime.now()
print("proc %d end processing at %s" %(rank,proc_end_time))
print("total processing time for proc %d is %s"%(rank,proc_end_time-proc_start_time))

# GATHER data : each proc send the result to rank 0
result_list = comm.gather(proc_result,root=0)

# rank 0 print the final result
if rank == 0:
    total_result = [0,Counter(),Counter()]
    listsToOne(result_list,total_result)
    printResult(total_result,query_word)

    end_time = datetime.datetime.now()
    print("end_time: ", end_time)
    IO_time = IO_end_time - start_time
    print("IO_time: ", IO_time)
    total_time = end_time - start_time
    print("total_time: ", total_time)

