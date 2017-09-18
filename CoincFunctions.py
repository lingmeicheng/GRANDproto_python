#!/home/sandra/anaconda3/bin/ipython

from SharedGlobals import *
import pickle
import numpy as np
import os.path
import struct #Interpret strings as packed binary data
from math import *

def EventTimeTableBuilder(nrun): #mettre Struct en argument au lieu de nrun quand termine

    #pour ouvrir struct au lieu de passer par DstBuilderLauncher.py
    runstr='R00'+str(nrun)
    setupfile="{0}setup_{1}".format(DST_DATA_PATH,runstr)
    with open(setupfile,'rb') as fs:
        unpick=pickle.Unpickler(fs)
        Struct=unpick.load()

    det=Struct['Setup']['Det'] #det is a list of ndets dictionnaries
    NbEventTotal=Struct['Setup']['TotalEvt']
    EventTimeTable=np.zeros((NbEventTotal,7))
    cpt_event=int(0)

    for i,elt in enumerate(det):        
        isScint=elt['isScint']
        if isScint:
            ib=ibuffs
        else:
            ib=ibuff

        for j in range(0,3):
            machine=Struct['Setup']['Det'][i]['Channels'][j]['Machine']  
            delay=Struct['Setup']['Det'][i]['Channels'][j]['Delay']
            machinestr=str(machine)
            print(machinestr)          

            filename="{0}{1}/{2}_A0{3}_time.bin".format(RAW_DATA_PATH,runstr,runstr,machinestr)
            if os.path.isfile(filename):
                with open(filename,'rb') as ft:
                    content=ft.read()
                    size=int(len(content)/4) #32bits data = 4 bytes data
                    content=struct.unpack('I'*size,content[0:len(content)]) #https://docs.python.org/2/library/struct.html
                    size=int(len(content)/ibufft) #4 times per machine event                    
                    if len(content)>=ibufft:
                        cpt_event_mach=int(1)  #not to give 0 to the first machine event id
                        while cpt_event_mach<=size:
                            
                            EventTimeTable[cpt_event,0]=(content[(cpt_event_mach-1)*ibufft+1]-1)*2**28 +(content[(cpt_event_mach-1)*ibufft+2]-1)*ib +content[(cpt_event_mach-1)*ibufft+3] -delay
                            EventTimeTable[cpt_event,1]=machine
                            EventTimeTable[cpt_event,2]=cpt_event_mach
                            EventTimeTable[cpt_event,3]=isScint
                            EventTimeTable[cpt_event,4]=content[(cpt_event_mach-1)*ibufft]
                            EventTimeTable[cpt_event,6]=elt['Id']

                            cpt_event_mach=cpt_event_mach+1
                            cpt_event=cpt_event+1
                    else:
                        print('File for machine '+machinestr+' is quite empty')

            else:
                print('Mei you file for machine '+machinestr)
            print(cpt_event)
                       
    abstimerefstart=EventTimeTable[:,0].min()*(1/FSAMPLING)
    abstimerefstop=EventTimeTable[:,0].max()*(1/FSAMPLING)
    totaltime=abstimerefstop-abstimerefstart
    nbloop=int(totaltime//TriggerTimeSpan)
    if totaltime%TriggerTimeSpan!=0:
        TrigRate=np.zeros((nbloop+1,len(det),3))
        TrigTime=np.zeros((nbloop+1,1))
        difftime=(totaltime%TriggerTimeSpan)/TriggerTimeSpan
    else:
        TrigRate=np.zeros((nbloop,len(det),3))
        TrigTime=np.zeros((nbloop,1))
        difftime=0

    order=np.argsort(EventTimeTable[:,0])
    EventTimeTable=EventTimeTable[order,:]
    time=EventTimeTable[:,0]*(1/FSAMPLING) #this is not an alias
    mach=EventTimeTable[:,1] #this is an alias
    #to copy an array, dont use a=b but a=b.copy(), otherwise it in only an alias

    '''loop=int(0)
    print('Trigger loop is running..')
    while loop<nbloop:
        if loop/100==floor(loop/100):
            print(loop)
        TrigTime[loop][0]=Struct['Setup']['RunTimeStart']+loop*TriggerTimeSpan
        for i,elt in enumerate(det):
            for j in range(0,3):
                machine=Struct['Setup']['Det'][i]['Channels'][j]['Machine']
                #nonzero gives the coordinates on the array that respect the condition,
                #with an array for the row coordinate and one for the column
                indmach=np.nonzero(mach==machine)
                indtimesup=np.nonzero(time>abstimerefstart+loop*TriggerTimeSpan)
                indtimeinf=np.nonzero(time<=abstimerefstart+(loop+1)*TriggerTimeSpan)
                indtime=np.intersect1d(indtimesup[0],indtimeinf[0])
                indtot=np.intersect1d(indtime,indmach[0])
                size=np.shape(indtot)[0] #row size of indtot
                TrigRate[loop][i][j]=size/TriggerTimeSpan
                EventTimeTable[indtot,5]=TrigRate[loop][i][j]
        loop=loop+1
    if difftime>0:
        TrigTime[nbloop][0]=Struct['Setup']['RunTimeStart']+nbloop*TriggerTimeSpan
        for i,elt in enumerate(det):
            for j in range(0,3):
                machine=Struct['Setup']['Det'][i]['Channels'][j]['Machine']
                #nonzero gives the coordinates on the array that respect the condition,
                #with an array for the row coordinate and one for the column
                indmach=np.nonzero(mach==machine)
                indtimesup=np.nonzero(time>abstimerefstart+nbloop*TriggerTimeSpan)
                indtimeinf=np.nonzero(time<=abstimerefstop)
                indtime=np.intersect1d(indtimesup[0],indtimeinf[0])
                indtot=np.intersect1d(indtime,indmach[0])
                size=np.shape(indtot)[0] #row size of indtot
                TrigRate[nbloop][i][j]=size/TriggerTimeSpan
                EventTimeTable[indtot,5]=TrigRate[nbloop][i][j]'''


    Struct['Setup']['InfosRun']['RawTrigRate']=TrigRate
    Struct['Setup']['InfosRun']['TrigTime']=TrigTime
    
    #saving the struct dictionnary in setupfile

    setupfile="{0}setup_{1}".format(DST_DATA_PATH,runstr)
    with open(setupfile,'wb') as fs:
        pick=pickle.Pickler(fs)
        pick.dump(Struct)

    print(np.shape(TrigRate))
    #print(TrigRate[1,:,:])
    print(EventTimeTable[0:15,0:7])

    return EventTimeTable, Struct
        



  





def CoincidenceFinder(nrun,EventTimeTable): #mettre Struct en argument au lieu de nrun quand termine

    #there is probably bugs in this function

    #pour ouvrir struct au lieu de passer par DstBuilderLauncher.py
    runstr='R00'+str(nrun)
    setupfile="{0}setup_{1}".format(DST_DATA_PATH,runstr)
    with open(setupfile,'rb') as fs:
        unpick=pickle.Unpickler(fs)
        Struct=unpick.load()

    TWinFactor=1.5
    duplication=0

    machines=[0,0,0] #for enumerate
    det=Struct['Setup']['Det'] #det is a list of ndets dictionnaries
    CoincNb=int(0)
    MultDet=[]
    MultSci=[]
    IDCoinc=[]
    IDDet=[]
    TagDet=[]
    UnixTime=[]
    EvtDet=[]
    StatusDet=[]
    TimeDet=[]
    TriggerRate=[]
    print(np.shape(EventTimeTable)[0])
    RefDet=int(0)
    while RefDet<np.shape(EventTimeTable)[0]: #number of rows

        NbDets=int(1)
        NbMachines=int(1)
        RefDetTime=EventTimeTable[RefDet,0]
        RefDetId=EventTimeTable[RefDet,6]
        machine=EventTimeTable[RefDet,1]
        CaracDetCoinc=np.zeros((len(det),7,3))
        CaracDetCoinc[0,0,:]=EventTimeTable[RefDet,6] #Id
        CaracDetCoinc[0,1,:]=EventTimeTable[RefDet,3] #isScint
        CaracDetCoinc[0,2,:]=EventTimeTable[RefDet,4] #unixtime

        for i,elt in enumerate(det):
            if elt['Id']==RefDetId:
                PosRefDet=np.array([elt['X'], elt['Y'], elt['Z']])
                for j,eltch in enumerate(machines):
                    if elt['Channels'][j]['Machine']==machine:
                        CaracDetCoinc[0,3,j]=EventTimeTable[RefDet,1] #machine
                        CaracDetCoinc[0,4,j]=EventTimeTable[RefDet,2] #machine event number
                        CaracDetCoinc[0,5,j]=EventTimeTable[RefDet,0] #triggertime in 5ns bins
                        CaracDetCoinc[0,6,j]=EventTimeTable[RefDet,5] #triggerate at the event moment



        TestDet=RefDet+1
        if TestDet+2*len(det)*3<np.shape(EventTimeTable)[0]: #number of rows
            cond=TestDet+2*len(det)*3
        if TestDet+2*len(det)*3>=np.shape(EventTimeTable)[0]: #number of rows
            cond=np.shape(EventTimeTable)[0]
        while TestDet<cond:
            TestDetTime=EventTimeTable[TestDet,0]
            TestDetId=EventTimeTable[TestDet,6]
            machine=EventTimeTable[TestDet,1]
            for i,elt in enumerate(det):
                if elt['Id']==TestDetId:
                    PosTestDet=np.array([elt['X'], elt['Y'], elt['Z']])
                    for j,eltch in enumerate(machines):
                        if elt['Channels'][j]['Machine']==machine:
                            indc=j
            DistDet=np.linalg.norm(PosRefDet-PosTestDet)
            DistMax=DistDet*FSAMPLING/C0 #number of 5ns bins necessary to light to cross DistDet
            if (TestDetTime-RefDetTime)<DistMax*TWinFactor+10:
                indd=np.nonzero(CaracDetCoinc[:,0,0]==TestDetId)[0]                
                if np.shape(indd)[0]>0: #if this detector has already been considered
                    if (np.shape(np.nonzero(CaracDetCoinc[indd[0],3,:]==machine)[0])[0])==0: #if the channel has not also already been considered
                        #print(machine, 'on entre ici parfois!', CaracDetCoinc[indd[0],3,:])
                        NbMachines=NbMachines+1 #increase of the number of machines for this coinc event                 
                        CaracDetCoinc[indd[0],3,indc]=EventTimeTable[TestDet,1] #machine
                        CaracDetCoinc[indd[0],4,indc]=EventTimeTable[TestDet,2] #machine event number
                        CaracDetCoinc[indd[0],5,indc]=EventTimeTable[TestDet,0] #triggertime in 5ns bins
                        CaracDetCoinc[indd[0],6,indc]=EventTimeTable[TestDet,5] #triggerate at the event moment
                        #print(CaracDetCoinc[indd[0],3,:])
                    else:
                        print('There is 1 duplication!!',TestDetTime,CaracDetCoinc[indd[0],5,indc])
                        duplication=duplication+1
                else:
                    NbMachines=NbMachines+1 #increase of the number of machines for this coinc event
                    NbDets=NbDets+1 #increase of the number of detectors for this coinc event
                    CaracDetCoinc[NbDets-1,0,:]=EventTimeTable[TestDet,6] #Id
                    CaracDetCoinc[NbDets-1,1,:]=EventTimeTable[TestDet,3] #isScint
                    CaracDetCoinc[NbDets-1,2,:]=EventTimeTable[TestDet,4] #unixtime
                    CaracDetCoinc[NbDets-1,3,indc]=EventTimeTable[TestDet,1] #machine
                    CaracDetCoinc[NbDets-1,4,indc]=EventTimeTable[TestDet,2] #machine event number
                    CaracDetCoinc[NbDets-1,5,indc]=EventTimeTable[TestDet,0] #triggertime in 5ns bins
                    CaracDetCoinc[NbDets-1,6,indc]=EventTimeTable[TestDet,5] #triggerate at the event moment
            else: #not anymore event machine in coinc for this coinc
                break
            TestDet=TestDet+1

        
        NbSci=np.shape(np.nonzero(CaracDetCoinc[:,1,1]==1)[0])[0]
        if AnalysisType==0:
            CoincSelection=NbDets>=5
        elif AnalysisType==1:
            #CoincSelection=(NbDets>=5 and NbSci>=2 and NbDets>NbSci)
            CoincSelection=(NbDets>=4 and NbSci>=1 and NbDets>NbSci)
        elif AnalysisType==2:
            CoincSelection=NbSci>=3

        if CoincSelection:
            CoincNb=CoincNb+1

            Carac2D=CaracDetCoinc[0:NbDets,0:3,0] #is an alias; Id,isScint,unixtime
            order=np.argsort(Carac2D[:,0])
            Carac2D=Carac2D[order,:]
            Carac3D=CaracDetCoinc[0:NbDets,3:7,0:3] #is an alias; machine, machine event number, triggertime in 5ns bins, triggerate at the event moment
            Carac3D=Carac3D[order,:,:]

            MultDet.append(NbDets)
            MultSci.append(NbSci)           
            IDCoinc.append(CoincNb)
            tag=[]
            ID=[]
            ut=[]
            evt=[]
            stat=[]
            time=[]
            trig=[]

            for i,elt in enumerate(det):
                evt.append([0,0,0])
                stat.append([0,0,0])
                time.append([0,0,0])
                trig.append([0,0,0])
                ind=np.nonzero(Carac2D[:,0]==elt['Id'])[0]
                if np.shape(ind)[0]>0:
                    tag.append(1)
                    ID.append(elt['Id'])
                    ut.append(Carac2D[ind,2][0])
                    for j,eltch in enumerate(machines):
                        evt[i][j]=int(Carac3D.tolist()[ind][1][j])
                        #print(evt[i][j])
                        stat[i][j]=Carac3D.tolist()[ind][0][j] #a arranger
                        time[i][j]=Carac3D.tolist()[ind][2][j]
                        trig[i][j]=Carac3D.tolist()[ind][3][j]              
                else:
                    tag.append(0)
                    ID.append(0)
                    ut.append(0)


            IDDet.append(ID)
            TagDet.append(tag)
            UnixTime.append(ut)
            EvtDet.append(evt)
            StatusDet.append(stat)
            TimeDet.append(time)
            TriggerRate.append(trig)
            #[[[1, 2], [2, 3], [3, 4]], [[4, 5], [5, 6], [6, 7]]]
            #2ligne 3 colonnes 2 deep
        
            RefDet=RefDet+NbMachines
        else:
            RefDet=RefDet+1
    print(CoincNb)
    print(duplication)

    #building dictionnaries hierarchy

    Coinc=dict()
    Det=dict()
    Coinc['Mult']=np.asarray(MultDet)
    Coinc['MultSci']=np.asarray(MultSci)
    Coinc['MultAnt']=Coinc['Mult']-Coinc['MultSci']
    Coinc['IdCoinc']=np.asarray(IDCoinc)
    Det['Tag']=np.asarray(TagDet)
    Det['UnixTime']=np.asarray(UnixTime)
    Det['Id']=np.asarray(IDDet)
    Det['Status']=np.asarray(StatusDet)
    Det['Evt']=np.asarray(EvtDet)
    Det['Time']=np.asarray(TimeDet)
    Det['TriggerRate']=np.asarray(TriggerRate)
    Coinc['Det']=Det
    Struct['Coinc']=Coinc
    Struct['Setup']['TotalCoinc']=CoincNb

    #saving the struct dictionnary in coincfile

    coincfile="{0}coinc_{1}".format(DST_DATA_PATH,runstr)
    with open(coincfile,'wb') as fs:
        pick=pickle.Pickler(fs)
        pick.dump(Struct)


    return Struct






if __name__ == "__main__": #si le module nest pas importe mais execute seul
    EventTimeTableBuilder(7005)





