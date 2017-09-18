#!/home/sandra/anaconda3/bin/ipython

from SharedGlobals import *
from ToolsFunctions import *
import pickle
import numpy as np
import os.path
import struct #Interpret strings as packed binary data
from math import *
#from scipy import stats
from scipy import interpolate
import matplotlib.pyplot as plt


def ConsecutiveCoincidenceFilter(nrun): #mettre Struct en argument au lieu de nrun quand termine

    #pour ouvrir struct au lieu de passer par DstBuilderLauncher.py
    runstr='R00'+str(nrun)
    coincfile="{0}coinc_{1}".format(DST_DATA_PATH,runstr)
    with open(coincfile,'rb') as fs:
        unpick=pickle.Unpickler(fs)
        Struct=unpick.load()

    NbCoinc=Struct['Setup']['TotalCoinc']
    Tag=Struct['Coinc']['Det']['Tag']
    Time=Struct['Coinc']['Det']['Time']
    TimeDiff=np.zeros(NbCoinc)
    EvtFlag=np.zeros(NbCoinc)
    Reject=np.zeros(NbCoinc)

    upperlimit=0.1 #sec

    cpt=0    
    while cpt<NbCoinc-1:
        Ref=cpt
        TagRef=np.nonzero(Tag[Ref,:])
        TagCons=np.nonzero(Tag[Ref+1,:])

        TimeRef=Time[Ref,TagRef[0],:]
        #print('TimeRef',TimeRef)
        TimeRef=min(TimeRef[TimeRef>0])
        #print('TimeRef',TimeRef)
        TimeCons=Time[Ref+1,TagCons[0],:]
        #print('TimeCons',TimeCons)
        TimeCons=min(TimeCons[TimeCons>0])
        #print('TimeCons',TimeCons)

        comAnt=np.intersect1d(TagRef[0],TagCons[0],True)
        TimeDiff[Ref]=TimeCons-TimeRef


        if len(comAnt)==0: #impossible for grandproto with 6antennas!
            EvtFlag[Ref]=1

        cpt=cpt+1

    TimeDiff=TimeDiff*TSAMPLING 

    if ConsCoincCompleteRejection:
        for i in range(0,NbCoinc-1):
            if TimeDiff[i]<upperlimit and EvtFlag[i]==0:
                Reject[i]=1
                Reject[i+1]=1
            #print(TimeDiff[i])
    else:
        print('Uncomplete consecutive coincidence rejection not implemented.')


    Rejection=dict()
    Rejection['ConsCoinc']=Reject
    Struct['Coinc']['Reject']=Rejection
    print(len(np.nonzero(Reject)[0])) #number of rejected events


    #saving the struct dictionnary in preprocessfile

    preprocessfile="{0}preprocess_{1}".format(DST_DATA_PATH,runstr)
    with open(preprocessfile,'wb') as fs:
        pick=pickle.Pickler(fs)
        pick.dump(Struct)


    return Struct









def RawFilter(nrun):
    
    #pour ouvrir struct au lieu de passer par DstBuilderLauncher.py
    runstr='R00'+str(nrun)
    preprocessfile="{0}preprocess_{1}".format(DST_DATA_PATH,runstr)
    with open(preprocessfile,'rb') as fs:
        unpick=pickle.Unpickler(fs)
        Struct=unpick.load()

    print(Struct['Setup']['TotalCoinc'])
    NbCoinc=Struct['Setup']['TotalCoinc']
    det=Struct['Setup']['Det'] #det is a list of ndets dictionnaries
    EvtId=Struct['Coinc']['Det']['Evt']

    settings=dict()
    settings['threshold']=threshold
    settings['granularity_box']=granularity_box
    settings['granularity_sampling']=granularity_sampling
    settings['max_out']=max_out
    settings['max_out_ToT']=max_out_ToT
    settings['max_block_ToT']=max_block_ToT
    settings['noise_mode']=noise_mode
    rawfilter=dict()
    rawfilter['settings']=settings
    Struct['rawfilter']=rawfilter

    MultRef=Struct['Coinc']['Mult']+0 #add a zero to avoid aliases
    Mu=np.zeros((NbCoinc,len(det),3))
    Sigma=np.zeros((NbCoinc,len(det),3))
    minraw=np.zeros((NbCoinc,len(det),3))-1
    maxraw=np.zeros((NbCoinc,len(det),3))-1
    sat=np.zeros((NbCoinc,len(det),3),int)
    nobox=np.zeros((NbCoinc,len(det),3))
    CaracEvt=np.zeros((NbCoinc,len(det),3),list)

    Granularity_box=int(round(settings['granularity_box']*FSAMPLING)) #3
    Granularity_sampling=int(round(settings['granularity_sampling']*FSAMPLING)) #8

    for i,elt in enumerate(det):
    #for i,elt in enumerate(det[0:1]):
        isScint=elt['isScint']
        if isScint:
            ib=ibuffs
        else:
            ib=ibuff

        tmu=np.arange(0,ib,1) #1023 values
        tmu=tmu*TSAMPLING*1e6 #tmu mightbe useless
        for j in range(0,3):
            machinestr=str(elt['Channels'][j]['Machine'])            
            filename="{0}{1}/{2}_A0{3}_data.bin".format(RAW_DATA_PATH,runstr,runstr,machinestr)
            print(elt['Id'],machinestr)
            if os.path.isfile(filename):
                ind=np.nonzero(EvtId[:,i,j])[0]
                print(len(ind))
                ngood=0
                with open(filename,'rb') as fd:
                    #content=fd.read()
                    #size=int(len(content)) #8bits data = 1 bytes data
                    #content=struct.unpack('B'*size,content) #https://docs.python.org/2/library/struct.html

                    for k in range(0,len(ind)):
                        thisEvt=EvtId[ind[k],i,j]-1
                        fd.seek(thisEvt*ib,0)
                        DataEvt=fd.read(ib)
                        DataEvt=struct.unpack('B'*ib,DataEvt)
                        DataEvt=np.asarray(DataEvt)
                        maxraw[ind[k],i,j]=max(DataEvt)
                        minraw[ind[k],i,j]=min(DataEvt)
                        if maxraw[ind[k],i,j]==255 or minraw[ind[k],i,j]==0:
                            sat[ind[k],i,j]=1

                        #if thisEvt+1==241 and i==5 and j==0:
                        #if thisEvt+1==353 and i==0 and j==2:
                        #    plt.plot(DataEvt)
                        #    plt.show()

                        KLOW=np.arange(0,round(ib/2-100),1,int) #0 to 411, 412 values
                        indmax=np.nonzero(max(DataEvt[KLOW])==DataEvt[KLOW])[0][0]
                        if max(DataEvt[KLOW])>155 and indmax>150:
                            KLOW=np.arange(0,indmax-10,1)

                        if settings['noise_mode']==1:
                            hist,edges=np.histogram(DataEvt[KLOW],bins=256,range=(0,256))# 0 values on the bin 0,...,255 values on the bin 255 (256th bin)
                            histmax=max(hist)
                            indmax=np.nonzero(hist==max(hist))[0][0]
                            Dx=np.unique(abs(edges[hist>0.5*histmax]-indmax))
                            if len(Dx<4):
                                mu=np.mean(DataEvt[KLOW])
                                sig=np.std(DataEvt[KLOW])
                            else:
                                n=round(Dx[-1]*2.5480) #3sigma
                                #n=round(Dx[-1]*2.5577) #variante
                                k=indmax+np.arange(-n,n+1,1)
                                pk=hist[k]/np.sum(hist[k])
                                mu=sum(k*pk)
                                sig=sqrt(sum((k-mu)**2*pk))
                        else:
                            mu=np.mean(DataEvt[KLOW])
                            sig=np.std(DataEvt[KLOW])

                        Mu[ind[k],i,j]=mu
                        Sigma[ind[k],i,j]=sig
                        

                        if isScint==0:
                            indWin=abs(DataEvt-mu)>=settings['threshold']*sig #true for values 5sigma above the mean, 0 else
                            Win=np.zeros(ib)
                            Win[indWin]=1
                            if thisEvt+1==353 and i==0 and j==2:
                                print(mu,sig)
                                print(Win[500:ib])
                            if np.sum(Win)>1:
                                Win=Agglomerate(Win, Granularity_sampling)
                            if thisEvt+1==353 and i==0 and j==2:
                                print(Win[500:ib])
                        else:
                            indWin=abs(DataEvt-mu)>=2*settings['threshold']*sig
                            Win=np.zeros(ib)
                            Win[indWin]=1
                            Win=Agglomerate(Win, Granularity_sampling)

                        S=np.nonzero(Win)[0]
                        if len(S)==0:
                            n_blocks=0
                            block_start=[]
                            block_end=[]
                        else:
                            dS=np.zeros(len(S))
                            dS[0]=2
                            dS[1:len(S)]=np.diff(S)
                            bg=np.nonzero(dS>1)[0]
                            ed=np.zeros(len(bg),int)
                            ed[0:len(bg)-1]=bg[1:len(bg)]-1
                            ed[-1]=len(S)-1
                            n_blocks=len(bg)
                            block_start=S[bg]
                            block_end=S[ed]
                            
                            for l in range(0,n_blocks):
                                ibef=block_start[l]-1
                                sel=np.nonzero(abs(DataEvt[max(1,ibef-Granularity_box):ibef+1]-mu)>=(settings['threshold']-2)*sig)[0] #3sigma
                                selind=ibef-Granularity_box+sel #works only if max is ibef-Granularity_box
                                selind=selind[selind>=0] #otherwise
                                Win[selind]=1

                                iaft=block_end[l]+1
                                sel=np.nonzero(abs(DataEvt[iaft:min(ib,iaft+Granularity_box+1)]-mu)>=(settings['threshold']-2)*sig)[0] #3sigma
                                selind=iaft+sel
                                selind=selind[selind<ib] #otherwise
                                Win[selind]=1
                        if thisEvt+1==353 and i==0 and j==2:
                            print(Win[500:ib])
                                
                        Win=Agglomerate(Win,Granularity_box)
                        if thisEvt+1==353 and i==0 and j==2:
                            print(Win[500:ib])
                        S=np.nonzero(Win)[0]
                        if len(S)==0:
                            n_blocks=0
                            block_start=[]
                            block_end=[]
                            block_len=[]
                            block_amp=[]
                            time_over_threshold=0
                            block_dt=[]
                        else:
                            dS=np.zeros(len(S))
                            dS[0]=2
                            dS[1:len(S)]=np.diff(S)
                            bg=np.nonzero(dS>1)[0]
                            ed=np.zeros(len(bg),int)
                            ed[0:len(bg)-1]=bg[1:len(bg)]-1
                            ed[-1]=len(S)-1
                            n_blocks=len(bg)
                            block_start=S[bg]
                            block_end=S[ed]

                            block_len=(block_end-block_start+1)*TSAMPLING
                            time_over_threshold=sum(block_len)
                            block_amp=np.zeros(len(block_len))
                            for l in range(0,n_blocks):
                                K=np.arange(block_start[l],block_end[l]+1,1,int)
                                block_amp[l]=max(abs(DataEvt[K]-mu))
                            if n_blocks>=2:
                                block_dt=np.zeros(n_blocks)
                                block_dt[0]=1
                                block_dt[1:len(block_dt)]=block_start[1:len(block_start)]-block_end[0:len(block_end)-1]
                            else:
                                block_dt=[1]
                        
                        if n_blocks>=1:
                            ngood=ngood+1
                            is_selected=1
                        else:
                            is_selected=0
                            MultRef[ind[k]]=MultRef[ind[k]]-1
                            nobox[ind[k],i,j]=1 #there is no boxes for this trace

                        cellEvt=[]
                        cellEvt.append(is_selected)
                        cellEvt.append(n_blocks)
                        cellEvt.append(time_over_threshold)
                        cellEvt.append(block_start)
                        cellEvt.append(block_end)
                        cellEvt.append(block_len)
                        cellEvt.append(block_amp)
                        cellEvt.append(block_dt)

                        CaracEvt[ind[k],i,j]=cellEvt

                    

    reject=np.zeros(len(MultRef))
    indr=np.nonzero(MultRef<4)[0]
    reject[indr]=1
    print(reject)

    Struct['Coinc']['Det']['Sigma']=Sigma
    Struct['Coinc']['Det']['Mu']=Mu
    Struct['Coinc']['Det']['MinRaw']=minraw
    Struct['Coinc']['Det']['MaxRaw']=maxraw
    Struct['Coinc']['Det']['Sat']=sat

    Struct['Coinc']['Reject']['RawFilter']=reject
    Struct['Coinc']['Reject']['NewMult']=MultRef
    Struct['Coinc']['Reject']['CaracEvt']=CaracEvt
    Struct['Coinc']['Reject']['NoBox']=nobox
    ind=np.nonzero(EvtId[:,1,1])[0]
    #print(Struct['Coinc']['Reject']['CaracEvt'][0,0,0])
    #print(Struct['Coinc']['Reject']['CaracEvt'][1,0,0])
    #print(Struct['Coinc']['Reject']['CaracEvt'][ind[0],1,1])
    #print(Struct['Coinc']['Reject']['CaracEvt'][ind[1],1,1])

    #saving the struct dictionnary in preprocessfile

    preprocessfile="{0}preprocess_{1}".format(DST_DATA_PATH,runstr)
    with open(preprocessfile,'wb') as fs:
        pick=pickle.Pickler(fs)
        pick.dump(Struct)

    return Struct









def CoincidenceFiltering(nrun):

    #pour ouvrir struct au lieu de passer par DstBuilderLauncher.py
    runstr='R00'+str(nrun)
    preprocessfile="{0}preprocess_{1}".format(DST_DATA_PATH,runstr)
    with open(preprocessfile,'rb') as fs:
        unpick=pickle.Unpickler(fs)
        preStruct=unpick.load()

    NbCoinc=preStruct['Setup']['TotalCoinc']
    CaracEvt=preStruct['Coinc']['Reject']['CaracEvt']
    EvtId=preStruct['Coinc']['Det']['Evt']+0 #add 0 to avoid aliases
    tag=preStruct['Coinc']['Det']['Tag']
    Sigma=preStruct['Coinc']['Det']['Sigma']
    det=preStruct['Setup']['Det'] #det is a list of ndets dictionnaries
    
    tagfilter=np.zeros(np.shape(EvtId))
    goodevts=np.zeros(NbCoinc,int)
    FiltResult=np.zeros(np.shape(EvtId))

    for i in range(0,NbCoinc):
        inddet=np.nonzero(tag[i,:])[0]
        for j in range(0,len(inddet)):            
            indch=np.nonzero(EvtId[i,inddet[j],:])[0]
            for k in range(0,len(indch)):
                filtresult=np.zeros((5,2))
                tagfilter[i,inddet[j],indch[k]]=1
                #print(CaracEvt[i,inddet[j],indch[k]])
                cellEvt=CaracEvt[i,inddet[j],indch[k]]
                #print(i,j,k,EvtId[i,inddet[j],indch[k]])
                #print(len(cellEvt),cellEvt)
                is_selected=cellEvt[0]
                n_blocks=cellEvt[1]
                #time_over_threshold=cellEvt[2]
                block_start=cellEvt[3]
                block_end=cellEvt[4]
                block_len=cellEvt[5]
                block_amp=cellEvt[6]
                block_dt=cellEvt[7]

                Threshold=threshold*Sigma[i,inddet[j],indch[k]]

                if n_blocks==0:
                    filtresult[4,0]=n_blocks
                    filtresult[4,1]=1
                    tagfilter[i,inddet[j],indch[k]]=0

                else:
                    block_mid=(block_start+block_end)/2
                    a=min(abs(block_mid-ibuff/2))
                    indcenter=np.nonzero(min(abs(block_mid-ibuff/2))==a)[0]

                    TotalPreTrig=(ibuff/2-block_start[indcenter])/FSAMPLING
                    if TotalPreTrig>max_pretrig_ToT:
                        filtresult[1,0]=TotalPreTrig
                        filtresult[1,1]=1

                    if block_len[indcenter]>max_block_ToT:
                        filtresult[2,0]=block_len[indcenter]
                        filtresult[2,1]=1
                        tagfilter[i,inddet[j],indch[k]]=0

                    if n_blocks>1:
                        time_over_threshold = sum(block_len)-block_len[indcenter]
                        if time_over_threshold>max_out_ToT:
                            filtresult[0,0]=time_over_threshold
                            filtresult[0,1]=1
                            tagfilter[i,inddet[j],indch[k]]=0

                    if n_blocks>1:
                        big=np.nonzero(block_amp>Threshold)[0]
                        close1=np.nonzero(block_start<ibuff/2)[0]
                        close2=np.nonzero(abs(block_mid-ibuff/2)<1e-6*FSAMPLING)[0]
                        large=np.nonzero(block_len>50e-9)[0]

                        close=np.intersect1d(close1,close2,True)
                        closebig=np.intersect1d(big,close,True)
                        bad=np.intersect1d(closebig,large,True)

                        ind=np.nonzero(bad==indcenter)[0]

                        if len(bad)-len(ind)>max_out:
                            filtresult[3,0]=len(bad)-len(ind)
                            filtresult[3,1]=1
                            tagfilter[i,inddet[j],indch[k]]=0

                if len(np.nonzero(filtresult[:,1]==1)[0]>0):
                    FiltResult[i,inddet[j],indch[k]]=np.nonzero(filtresult[:,1]==1)[0][0]+1 #reccord first cause of rejection
                if preStruct['Coinc']['IdCoinc'][i]==12:
                    print(cellEvt,indcenter)



        tagfinal=tagfilter[i,:,:].sum(1)
        indmult=np.nonzero(tagfinal>0)[0]
        nomult=np.nonzero(tagfinal==0)[0]
        noch=np.nonzero(tagfilter[i,:,:]==0)  
        #print(noch)

        if preStruct['Coinc']['IdCoinc'][i]==12:
            print(FiltResult[i,:,:])
            print(tagfilter[i,:,:])
            print(EvtId[i,:,:])
        if len(indmult)>3:
            goodevts[i]=1

            preStruct['Coinc']['Mult'][i]=len(indmult)
            preStruct['Coinc']['MultAnt'][i]=len(indmult)
            preStruct['Coinc']['Det']['Id'][i,nomult]=0
            preStruct['Coinc']['Det']['Tag'][i,nomult]=0
            preStruct['Coinc']['Det']['UnixTime'][i,nomult]=0
            preStruct['Coinc']['Det']['Status'][i,noch[0],noch[1]]=preStruct['Coinc']['Det']['Status'][i,noch[0],noch[1]]+2
            preStruct['Coinc']['Det']['Evt'][i,noch[0],noch[1]]=0
            preStruct['Coinc']['Det']['Time'][i,noch[0],noch[1]]=0 
            preStruct['Coinc']['Det']['TriggerRate'][i,noch[0],noch[1]]=0


            preStruct['Coinc']['Det']['Sigma'][i,noch[0],noch[1]]=0
            preStruct['Coinc']['Det']['Mu'][i,noch[0],noch[1]]=0
            preStruct['Coinc']['Det']['MinRaw'][i,noch[0],noch[1]]=-1
            preStruct['Coinc']['Det']['MaxRaw'][i,noch[0],noch[1]]=-1
            preStruct['Coinc']['Det']['Sat'][i,noch[0],noch[1]]=0

    preStruct['Coinc']['Det']['FiltResult']=FiltResult

    indgood=np.nonzero(goodevts)[0]
    print(indgood)
    for i,elt in enumerate(det):            
        for j in range(0,3):
            machinestr=str(elt['Channels'][j]['Machine']) 
            ind=np.nonzero(EvtId[:,i,j])[0] #EvtId is not an alias of th correspondng prestruct field
            if len(ind)!=0:
                nok=len(np.nonzero(FiltResult[ind,i,j]==0)[0])
                nout=len(np.nonzero(FiltResult[ind,i,j]==1)[0])
                npre=len(np.nonzero(FiltResult[ind,i,j]==2)[0])
                ncen=len(np.nonzero(FiltResult[ind,i,j]==3)[0])    
                nrep=len(np.nonzero(FiltResult[ind,i,j]==4)[0])
                nnob=len(np.nonzero(FiltResult[ind,i,j]==5)[0]) 
                nno=nout+npre+ncen+nrep+nnob
                ntot=len(ind);

                print('Machine '+machinestr, 'in',ntot,'coincs:',nok,'events OK (',nok/ntot*100,' pc)',nno, 'rejected (',nno/ntot*100,' pc)')
                print(nno, 'events rejected:',ncen,'too long central;',npre,'too early central;',nout,'too long outside;',nrep,'pulses outside;',nnob,'no box')
                print(len(np.nonzero(preStruct['Coinc']['Det']['Evt'][indgood,i,j])[0]))


    #new struct only with good coincs and good traces

    Struct=dict()
    Struct['Setup']=preStruct['Setup']
    Struct['Setup']['TotalCoinc']=len(indgood)

    Coinc=dict()
    Struct['Coinc']=Coinc

    Struct['Coinc']['Mult']=preStruct['Coinc']['Mult'][indgood]
    Struct['Coinc']['MultAnt']=preStruct['Coinc']['MultAnt'][indgood]
    Struct['Coinc']['MultSci']=preStruct['Coinc']['MultSci'][indgood]
    Struct['Coinc']['IdCoinc']=preStruct['Coinc']['IdCoinc'][indgood]

    Det=dict()
    Struct['Coinc']['Det']=Det

    Struct['Coinc']['Det']['Id']=preStruct['Coinc']['Det']['Id'][indgood,:]
    Struct['Coinc']['Det']['Tag']=preStruct['Coinc']['Det']['Tag'][indgood,:]
    Struct['Coinc']['Det']['UnixTime']=preStruct['Coinc']['Det']['UnixTime'][indgood,:]
    Struct['Coinc']['Det']['Status']=preStruct['Coinc']['Det']['Status'][indgood,:,:]
    Struct['Coinc']['Det']['Evt']=preStruct['Coinc']['Det']['Evt'][indgood,:,:]
    Struct['Coinc']['Det']['Time']=preStruct['Coinc']['Det']['Time'][indgood,:,:]
    Struct['Coinc']['Det']['TriggerRate']=preStruct['Coinc']['Det']['TriggerRate'][indgood,:,:]

    Struct['Coinc']['Det']['FiltResult']=preStruct['Coinc']['Det']['FiltResult'][indgood,:,:]

    Struct['Coinc']['Det']['Sigma']=preStruct['Coinc']['Det']['Sigma'][indgood,:,:]
    Struct['Coinc']['Det']['Mu']=preStruct['Coinc']['Det']['Mu'][indgood,:,:]
    Struct['Coinc']['Det']['MinRaw']=preStruct['Coinc']['Det']['MinRaw'][indgood,:,:]
    Struct['Coinc']['Det']['MaxRaw']=preStruct['Coinc']['Det']['MaxRaw'][indgood,:,:]
    Struct['Coinc']['Det']['Sat']=preStruct['Coinc']['Det']['Sat'][indgood,:,:]

    print(len(np.nonzero(Struct['Coinc']['Det']['Evt'])[0]))

    #saving the struct dictionnary in dstfile

    dstfile="{0}dst_{1}".format(DST_DATA_PATH,runstr)
    with open(dstfile,'wb') as fs:
        pick=pickle.Pickler(fs)
        pick.dump(Struct)

    return Struct










def TriggerTimeBuilder(nrun):  #mettre Struct en argument au lieu de nrun quand termine

    #pour ouvrir struct au lieu de passer par DstBuilderLauncher.py
    runstr='R00'+str(nrun)
    dstfile="{0}dst_{1}".format(DST_DATA_PATH,runstr)
    with open(dstfile,'rb') as fs:
        unpick=pickle.Unpickler(fs)
        Struct=unpick.load()

    print(Struct['Setup']['TotalCoinc'],len(np.nonzero(Struct['Coinc']['Det']['Evt'])[0]))
    det=Struct['Setup']['Det'] #det is a list of ndets dictionnaries
    EvtId=Struct['Coinc']['Det']['Evt']
    time=Struct['Coinc']['Det']['Time'] #time of the time data file (in 5ns bins)
    AmpMax=np.zeros(np.shape(time))-1
    TrigTime=np.zeros(np.shape(time))-1
    NbCoinc=Struct['Setup']['TotalCoinc']
    fover=10
    deltat=10 #real bins
    deltac=fover*deltat #100  

    for i,elt in enumerate(det): 
        isScint=elt['isScint']
        if isScint:
            ib=ibuffs
        else:
            ib=ibuff
        trig=int(ib/2)-1 #511 (512th bin)
        trigc=trig*fover #5110
        ts=np.arange(0,ib,1) #1023 values
        ts=ts*TSAMPLING
        tsover=np.arange(0,ib-1+1./fover,1./fover) #10231 values
        tsover=tsover*TSAMPLING
        for j in range(0,3):
            machinestr=str(elt['Channels'][j]['Machine'])            
            filename="{0}{1}/{2}_A0{3}_data.bin".format(RAW_DATA_PATH,runstr,runstr,machinestr)
            print(elt['Id'],machinestr)
            if os.path.isfile(filename):
                ind=np.nonzero(EvtId[:,i,j])[0]
                print(len(ind))
                with open(filename,'rb') as fd:
                    #content=fd.read()
                    #size=int(len(content)) #8bits data = 1 bytes data
                    #content=struct.unpack('B'*size,content) #https://docs.python.org/2/library/struct.html                
                    for k in range(0,len(ind)):
                        thisEvt=EvtId[ind[k],i,j]-1
                        fd.seek(thisEvt*ib,0)
                        DataEvt=fd.read(ib)
                        DataEvt=struct.unpack('B'*ib,DataEvt)
                        DataEvt=np.asarray(DataEvt)
                        DataEvt=DataEvt*SCALE
                        v=PassBand(DataEvt,ts,FREQMIN,FREQMAX)
                        tck=interpolate.splrep(ts,v)
                        vcor=interpolate.splev(tsover,tck)
                        #plt.plot(tsover,vcor,ts,v)
                        #plt.show()
                        v2=vcor[max(0,trigc-deltac):min(trigc+deltac,len(vcor))] #5010 to 5210 (200 values) -> 501 to 521
                        vpeakpeak=max(v2)-min(v2)
                        indmax=np.nonzero(v2==max(v2))[0][0]+trigc-deltac
                        indmin=np.nonzero(v2==min(v2))[0][0]+trigc-deltac
                        tmoy=(indmax+indmin)/2
                        AmpMax[ind[k],i,j]=vpeakpeak
                        TrigTime[ind[k],i,j]=time[ind[k],i,j]-trig+tmoy/fover
                        #print(TrigTime[ind[k],i,j])

    for i in range(0,NbCoinc):
        ind=np.nonzero(EvtId[i,:,:])
        TrigTime[i,ind[0][:],ind[1][:]]=TrigTime[i,ind[0][:],ind[1][:]]-min(TrigTime[i,ind[0][:],ind[1][:]])
    #print(TrigTime[NbCoinc-1,ind[0][:],ind[1][:]])

        if max(TrigTime[i,ind[0][:],ind[1][:]])>1000:
            print('TrigTime calculation error.')





    Struct['Coinc']['Det']['AmpMax']=AmpMax
    Struct['Coinc']['Det']['TrigTime']=TrigTime


    #saving the struct dictionnary in dstfile

    dstfile="{0}dst_{1}".format(DST_DATA_PATH,runstr)
    with open(dstfile,'wb') as fs:
        pick=pickle.Pickler(fs)
        pick.dump(Struct)

    return Struct





















if __name__ == "__main__": #si le module nest pas importe mais execute seul
    RawFilter(7005)
    #TriggerTimeBuilder(7005)
    #CoincidenceFiltering(7005)

