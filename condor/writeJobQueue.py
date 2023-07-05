import os,sys

indir = sys.argv[1]
files = []
for f in os.listdir(sys.argv[1]):
    if "output" in f and ".root" in f and os.path.isfile("%s/%s"%(sys.argv[1],f)):
        files.append("%s/%s"%(sys.argv[1],f))

fout = open("queue.txt","w")
for i in range(len(files)):
    fout.write("$ENV(SCOUTINGINPUTDIR) $ENV(SCOUTINGOUTPUTDIR) %d 1\n"%i)
