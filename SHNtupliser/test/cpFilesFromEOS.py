
import argparse

parser = argparse.ArgumentParser(description='gives the files and optional picks the events given')
parser.add_argument('--dataset',help='name of dataset')
parser.add_argument('--runnr',help='runnr')
parser.add_argument('--lumis',help='lumis seperated by :')
parser.add_argument('--outputDir',default="./",help='output directory')
args = parser.parse_args()
print args

from SHarper.SHNtupliser.dasFileQuery import *

eventList=[]
Event.dataset = args.dataset
import re
commentRE = re.compile (r'#.+$')
for lumi in args.lumis.split(":"):
    eventStr=str(args.runnr)+":"+lumi+":1"
    try:
        event = Event (eventStr)
    except:
        print "Skipping '%s'." % inputCleaned
    eventList.append(event)

filesSet=set([])
fileCache={}
for event in eventList:
    for filename in getFileNamesCache (event,fileCache):
        filesSet.add(filename)

for filename in filesSet:
    fileLumis = sorted(getLumis(filename)[0])
    
    fileData={'dataset' : args.dataset.split("/")[1],'runnr' : args.runnr,'minLS' : fileLumis[0],'maxLS' : fileLumis[-1]}
    outputFilename = "{d[dataset]}_{d[runnr]}_LS{d[minLS]}to{d[maxLS]}.root".format(d=fileData)
    print fileLumis,filename,outputFilename
    
    fullEOSFilename="/eos/cms/tier0/"+filename
    
    import subprocess
    out,err = subprocess.Popen(['/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select','cp',fullEOSFilename,args.outputDir+"/"+outputFilename],stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()
    for line in err.splitline():
        print line

