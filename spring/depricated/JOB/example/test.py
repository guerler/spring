#!/opt/anaconda-2.0/bin/python
from time import sleep
import os

def mkdir(directory):
    #make directory
    try:
        os.mkdir(directory)
    except:
        pass


username = 'bgovi'
fout_path = '/nfs/amino-home/'+username+'/test/'

#mkdirs
mkdir(fout_path)
mkdir("/tmp/"+username+"/")
mkdir("/tmp/"+username+"/example")

#write to file job_test a simp
os.chdir("/tmp/"+username+"/example")
fout = open(fout_path+"job_test.txt",'w')
for i in xrange(0,10):
    print i
    fout.write(str(i) +'\n')
    sleep(1)

sleep(1)
fout.close()
