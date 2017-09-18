#!/usr/bin/env python
#!/home/sandra/anaconda3/bin/ipython

#===============================================================================
# Grid Engine steering options
#===============================================================================
## Set the name
#$ -N dstbuilder

## Submit job under TREND group
#$ -P P_trend

## Merge the stdout et stderr in a single file
#$ -j y

## Files .e et .o copied to current working directory
#$ -cwd

## Notify stop and kill signals before issuing them.
#$ -notify

## CPU time
#$ -l ct=10:00:00

## Memory
#$ -l vmem=20.0G

## Disk space
#$ -l fsize=1.0G

## Request /sps/:
## -l sps=1

#===============================================================================




import subprocess #to launch linux commands
import os
import sys


IRODS_PATH='/trend/home/trirods/data/raw/'


from SharedGlobals import *
from SetupFunctions import *
from CoincFunctions import *
from ReconsFunctions import *
from ProcessingFunctions import *


def DstBuilderLauncher(*args):
    result=[]
    for arg in args:
        print("arguments:", arg)
        runstr='R00'+str(arg)
        dirname="{0}{1}/".format(RAW_DATA_PATH,runstr)
        print(dirname)
        if os.path.isdir(dirname):
            result.append(DstBuilder(arg))
    return result


def DstBuilder(nrun):
    runstr='R00'+str(nrun)

    print('RunSetupBuilder...')
    Struct=RunSetupBuilder(nrun) #we can erase all the analysis type arguments because it is in sharedglobals
    if Struct['Setup']['TotalEvt']==0: #if there is no event
        return Struct

    print('EventTimeTableBuilder...')
    EventTimeTable,Struct=EventTimeTableBuilder(nrun)
    print('CoincidenceFinder...')
    Struct=CoincidenceFinder(nrun,EventTimeTable)    
    if Struct['Setup']['TotalCoinc']==0: #if there is no event
        return Struct

    if 1:
        print('ConsecutiveCoincidenceFilter...')
        Struct=ConsecutiveCoincidenceFilter(nrun)
        if QuickReject==1:
            goodevts=np.nonzero(Struct['Coinc']['Reject']['ConsCoinc']==0)[0]
            Struct['Setup']['TotalCoinc']=len(goodevts)
            print(len(goodevts))
            Struct['Coinc']['Mult']=Struct['Coinc']['Mult'][goodevts]
            Struct['Coinc']['MultAnt']=Struct['Coinc']['MultAnt'][goodevts]
            Struct['Coinc']['MultSci']=Struct['Coinc']['MultSci'][goodevts]
            Struct['Coinc']['IdCoinc']=Struct['Coinc']['IdCoinc'][goodevts]
            Struct['Coinc']['Det']['Id']=Struct['Coinc']['Det']['Id'][goodevts,:]
            Struct['Coinc']['Det']['Tag']=Struct['Coinc']['Det']['Tag'][goodevts,:]
            Struct['Coinc']['Det']['UnixTime']=Struct['Coinc']['Det']['UnixTime'][goodevts,:]
            Struct['Coinc']['Det']['Status']=Struct['Coinc']['Det']['Status'][goodevts,:,:]
            Struct['Coinc']['Det']['Evt']=Struct['Coinc']['Det']['Evt'][goodevts,:,:]
            Struct['Coinc']['Det']['Time']=Struct['Coinc']['Det']['Time'][goodevts,:,:]
            Struct['Coinc']['Det']['TriggerRate']=Struct['Coinc']['Det']['TriggerRate'][goodevts,:,:]

    
    preprocessfile="{0}preprocess_{1}".format(DST_DATA_PATH,runstr)
    with open(preprocessfile,'wb') as fs:
        pick=pickle.Pickler(fs)
        pick.dump(Struct)

    print('RawFilter...')
    Struct=RawFilter(nrun)
    if AnalysisType==0:
        print('CoincidenceFiltering...')
        Struct=CoincidenceFiltering(nrun)
    else:
        dstfile="{0}dst_{1}".format(DST_DATA_PATH,runstr)
        with open(dstfile,'wb') as fs:
            pick=pickle.Pickler(fs)
            pick.dump(Struct)

    print('TriggerTimeBuilder...')
    Sruct=TriggerTimeBuilder(nrun)
    print('CoinctableFileBuilder...')
    Struct=CoinctableFileBuilder(nrun)
    if AnalysisType==0:
        reconsEx = CXX_PATH+'recons'
    elif AnalysisType==1:
        reconsEx = CXX_PATH+'reconsHybrid'
    elif AnalysisType==2:
        reconsEx = CXX_PATH+'reconsScint'
    print('Recons...')
    p=subprocess.Popen('export TREND_TEXT_PATH='+TEXT_PATH+'; '+reconsEx+' '+str(nrun)+' '+str(CORREL), shell=True)
    p.wait()

    print('ReconsLoader...')
    Struct=ReconsLoader(nrun)
    print('And finally they lived happy. The end!')
    return Struct



if __name__ == "__main__": #si le module nest pas importe mais execute seul
    print(len(sys.argv))
    for i in range(0,len(sys.argv)-1):
        pdir=subprocess.Popen('iget -vfr '+IRODS_PATH+'R00'+sys.argv[i+1]+' '+RAW_DATA_PATH, shell=True) 
        pdir.wait()
        DstBuilderLauncher(sys.argv[i+1]) #we can input here as many arguments as needed





