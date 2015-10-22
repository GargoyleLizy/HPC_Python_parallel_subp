import threading

def show(one_str):
    for onechar in one_str:
        print onechar



server_ip = ['123','345','678']
work = 'echo'

thds=[]
for ip in server_ip:
    thd= threading.Thread(target=show,args=(ip,)  )
    thds.append(thd)

for thd in thds:
    thd.start()

for thd in thds:
    thd.join()


