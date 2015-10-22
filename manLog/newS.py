import SocketServer
import socket,threading
import sys,os
import argparse

from ..Config.taskmaster import taskmaster
#from ..Config.cloudconfig import def_config
from ..Config.cloudconfig import *

# Initialization the configuration information

ab_cloud_dir=None


class EchoRequestHandler(SocketServer.BaseRequestHandler):
    data = 1
    CONFIG = def_config()
    SSTAT = server_status()
    PROTOC = protocol()

    def setup(self,temp_cloud_dir):
        # ****** A lot of checking/Preparation here *****
        # check if the project directories are set correctly
        #temp_cloud_dir = args.cdir
        (ab_cloud_dir,ab_projs_dir,projs_list) = self.check_projs_installed(temp_cloud_dir)
        if(ab_projs_dir == None ):
            print("project directory on the server side is not available. quit")
            exit()
        # SSAT record info
        self.SSTAT.set_ab_cloud_dir(ab_cloud_dir)
        self.SSTAT.set_ab_projs_dir(ab_projs_dir)
        self.SSTAT.set_projs_list(projs_list)

        # check the log of taks done 
        last_log_idx = self.SSTAT.check_taskslog()
        print('last job index is %d'%last_log_idx)
        return SocketServer.BaseRequestHandler.setup(self)

    # ***** some helper functions.*******
    # check if the server side projects are set appropriately. 
    def check_projs_installed(self,temp_cloud_dir):
        ab_cloud_dir = os.path.abspath(temp_cloud_dir)
        os.chdir(ab_cloud_dir)
        # check if Projs exist
        top_dirs=[ name for name in os.listdir('./') if os.path.isdir(os.path.join('./', name))]
        print(ab_cloud_dir)
        print("Top modules: %s "% top_dirs)
        if self.SSTAT.projs_dir in top_dirs:
            print('Find Projects directory. List available projects below:')
            ab_projs_dir = ab_cloud_dir + '/' + self.SSTAT.projs_dir
            os.chdir(ab_projs_dir)

            top_projs=[ name for name in os.listdir('./') 
                    if os.path.isdir(os.path.join('./', name))]
            print('Avail Projs: %s'%top_projs)
            #return (None,None,None)
            return (ab_cloud_dir,ab_projs_dir,top_projs) 
        else:
            return (None,None,None)
    



    def handle(self):
        # Echo the back to the client
        data = self.request.recv(1024)
        self.data+=1
        self.request.send(self.SSTAT.get_curent_status)
        return

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-cd','--cdir')
    parser.add_argument('-t','--test')
    args=parser.parse_args()


    address = ('localhost', 0) # let the kernel give us a port
    server = SocketServer.TCPServer(address, EchoRequestHandler )
    ip, port = server.server_address # find out what port we were given

    t = threading.Thread(target=server.serve_forever)
    t.setDaemon(True) # don't hang on exit
    t.start()

    # Connect to the server
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.connect((ip, port))

    # Send the data
    message = 'Hello, world'
    print 'Sending : "%s"' % message
    len_sent = s.send(message)

    # Receive a response
    response = s.recv(len_sent)
    s.send(message)
    response = s.recv(1024)
    print 'Received: "%s"' % response

    # Clean up
    s.close()
    server.socket.close()
