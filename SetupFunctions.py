#!/home/sandra/anaconda3/bin/ipython

from SharedGlobals import *
import struct #Interpret strings as packed binary data
import os.path
import numpy as np
import pickle
from math import *
 

def RunSetupBuilder(nrun):

    runstr='R00'+str(nrun)

    if AnalysisType==0:
        dets=ALL3DANTS
    elif AnalysisType==1:
        dets=ALL3DDETS
    elif AnalysisType==2:
        dets=ALL3DSCINTS

    posDet, podN, machines, delay, cable, isSci=SetupCharacteristics_GRANDproto(nrun,dets)

    #statement of main dictionnaries and variables

    Struct=dict()
    Setup=dict()
    Det=[]
    for elt in dets:
        Det.append(0)
    for i,elt in enumerate(dets):
        Det[i]=dict() #Det is a list of ndets dictionnaries
    InfosRun=dict()
    TimeStart=np.zeros((len(dets),3))
    TimeStop=np.zeros((len(dets),3))
    TotalEvt=int(0)
    RunTimeStart=0
    RunTimeStop=0

    #filling dictionnaries with SetupCharacteristics_GRANDproto results and time files

    for i,elt in enumerate(dets):
        Det[i]['Id']=elt
        Det[i]['isScint']=isSci[elt]
        Det[i]['X']=posDet[elt][0]
        Det[i]['Y']=posDet[elt][1]
        Det[i]['Z']=posDet[elt][2]

        Channels=[0,0,0] #if 3 channels
        for j,eltch in enumerate(Channels):
            Channels[j]=dict() #Channels is a list of 3 dictionnaries
            Channels[j]['Machine']=machines[elt][j]
            Channels[j]['PodNb']=podN[elt][j]
            Channels[j]['Cable']=cable[elt][j]
            Channels[j]['Delay']=delay[elt][j]

            #reading the time file of each machine
            
            machinestr=str(machines[elt][j])            
            filename="{0}{1}/{2}_A0{3}_time.bin".format(RAW_DATA_PATH,runstr,runstr,machinestr)

            if os.path.isfile(filename):
                with open(filename,'rb') as ft:
                    content=ft.read()
                    size=int(len(content)/4) #32bits data = 4 bytes data
                    content=struct.unpack('I'*size,content[0:len(content)]) #https://docs.python.org/2/library/struct.html
                    size=int(len(content)/ibufft) #4 times per machine event
                    TotalEvt=TotalEvt+size # raise of the total number of machine events
                    if len(content)>=ibufft:
                        Channels[j]['TimeStart']=content[0]
                        Channels[j]['TimeStop']=content[(size-1)*ibufft]
                        TimeStart[i][j]=content[0]
                        TimeStop[i][j]=content[(size-1)*ibufft]
                    else:
                        print('File for machine '+machinestr+' is quite empty')
                        Channels[j]['TimeStart']=0
                        Channels[j]['TimeStop']=0
            else:
                print('Mei you file for machine '+machinestr)
                Channels[j]['TimeStart']=0
                Channels[j]['TimeStop']=0

        Det[i]['Channels']=Channels

    #computation of RunTimeStart and RunTimeStop
    if len(np.nonzero(TimeStart)[0]) !=0:
        ts=TimeStart.flatten(1)
        ts=np.sort(ts)
        year=1970+(ts/(365.25*24*3600)) #rough
        ts=ts[year>=2015]
        diff=np.diff(ts)
        add=np.array([0])
        diff=np.hstack((diff,add))
        ts=ts[diff<60]
        RunTimeStart=ts.min()
        
        ts=TimeStop.flatten(1)
        ts=np.sort(ts)
        year=1970+(ts/(365.25*24*3600)) #rough
        ts=ts[year>=2015]
        diff=np.diff(ts)
        diff=np.hstack((diff,add))
        ts=ts[diff<60]
        RunTimeStop=ts.max()

    #etablishment of dictionnaries hierarchy

    InfosRun['TimeStart']=TimeStart #array
    InfosRun['TimeStop']=TimeStop #array
    Setup['Det']=Det #list of dictionnaries
    Setup['InfosRun']=InfosRun #dictionnary
    Setup['Run']=nrun #scalar
    Setup['TotalEvt']=TotalEvt #scalar
    Setup['RunTimeStart']=RunTimeStart #scalar
    Setup['RunTimeStop']=RunTimeStop #scalar
    Struct['Setup']=Setup

    #print(Struct['Setup']['Det'][1]['Channels'][0]['TimeStart'])
    #print(content[len(content)-100:len(content)])

    #saving the struct dictionnary in setupfile

    setupfile="{0}setup_{1}".format(DST_DATA_PATH,runstr)
    with open(setupfile,'wb') as fs:
        pick=pickle.Pickler(fs)
        pick.dump(Struct)

    return Struct



def SetupCharacteristics_GRANDproto(nrun,dets):

    if nrun<6700:
        print('SetupCharacteristics_GRANDproto not valid for runs below R6700')
        return
    
    posDet=dict()
    podN=dict()
    machines=dict()
    cable=dict()
    fiberdelay=dict()
    isSci=dict()
    delay=dict()

    for i in dets:
        ## Antennas 
        if i==1:          # Former B4
            posDet[i]=[19.0,    441.9,    2634.0] 
            podN[i]= [151, 151, 151]
            machines[i]= [129, 131, 132]
            if nrun>7149:
                machines[i]= [128, 131, 132]
            cable[i]= [+14, +14, +12]  # 9.5+9.2 for X, same for Y, 9.5+6.2 for Z
            fiberdelay[i]= [-70, -71 , -71] # 
            isSci[i]= 0
            delay[i]=[0,0,0]
            for j, elt in enumerate(delay[i]):
                delay[i][j]= cable[i][j]+fiberdelay[i][j]

        elif i==2: # C5
            posDet[i]=[152.8,466.3,2633.6]
            #pos_ant=[147, 461, 2640]  Prelim meas.
            podN[i]= [150, 150, 150]
            machines[i]= [134,135,136]
            cable[i]= [122, 124, 121] #R6899, 6902, 6903
            fiberdelay[i]= [-99, -99 , -99] #R6886, 6888, 6889  
            isSci[i]= 0
            delay[i]=[0,0,0]
            j=int(0)
            for j, elt in enumerate(delay[i]):
                delay[i][j]= cable[i][j]+fiberdelay[i][j]

        elif i==3:  # Former C3 
            posDet[i]=[-66.6,    477.7,    2636.1]
            podN[i]= [149,149,149]
            machines[i]= [137,138,140]
            cable[i]= [+56,+61,+56] #R6906, 6907, 6908
            fiberdelay[i]= [-122, -124 , -122] #R6905,6910,6909     
            isSci[i]= 0
            delay[i]=[0,0,0]
            j=int(0)
            for j, elt in enumerate(delay[i]):
                delay[i][j]= cable[i][j]+fiberdelay[i][j]

        elif i==4:  # Former C4    
            posDet[i]=[-6.7,    512.2,    2636.2]
            #pos_ant=[-10, 518, 2639]  # Prelim meas.
            podN[i]= [148,148,148]
            machines[i]= [151,154,155]
            cable[i]= [+8,+8,+8] #Mesure directe: cable[i]= 10m
            fiberdelay[i]= [-162, -160, -162] #R6911-6913                     
            isSci[i]= 0
            delay[i]=[0,0,0]
            j=int(0)
            for j, elt in enumerate(delay[i]):
                delay[i][j]= cable[i][j]+fiberdelay[i][j]

        elif i==5:  # Former A5 - 2D antenna
            posDet[i]=[-42.0,    635.0,    2640.1]
            podN[i]= [147, 146, 0]
            machines[i]= [148, 149, 0]
            cable[i]= [44.5, 81, 0]  # R6974, R6967
            fiberdelay[i]= [-298, -347.5 , 0] # R6922, R6938
            isSci[i]= 0
            delay[i]=[0,0,0]
            j=int(0)
            for j, elt in enumerate(delay[i]):
                delay[i][j]= cable[i][j]+fiberdelay[i][j]

        elif i==6:  # Former A3
            posDet[i]=[50.9,    657.6,    2640.3]
            podN[i]= [146, 146, 146]
            machines[i]= [156, 157, 158]
            cable[i]= [72.5, 70, 69.5]  #R6971,R6969,R6973 
            fiberdelay[i]= [-347.5, -348.5 , -347.5] #R6935-7
            isSci[i]= 0
            delay[i]=[0,0,0]
            j=int(0)
            for j, elt in enumerate(delay[i]):
                delay[i][j]= cable[i][j]+fiberdelay[i][j]
            
        ## Scints
        elif i==111:# S1 (West)
            posDet[i]=[-74.7,    391.2,    2633.7]
            podN[i]= [153,0,0]
            machines[i]= [111,0,0]
            cable[i]= [68,0,0] #R6871
            fiberdelay[i]= [1,0,0] #R6867
            isSci[i]= 1
            delay[i]=[0,0,0]
            for j, elt in enumerate(delay[i]):
                delay[i][j]= cable[i][j]+fiberdelay[i][j]

        elif i==112:# S2 (Center)
            posDet[i]=[15.4,    395.9,    2632.9]
            podN[i]= [153,0,0]
            machines[i]= [112,0,0]
            cable[i]= [+32,0,0 ] #R6872
            fiberdelay[i]= [0,0,0] #R6869   
            isSci[i]= 1
            delay[i]=[0,0,0]
            j=int(0)
            for j, elt in enumerate(delay[i]):
                delay[i][j]= cable[i][j]+fiberdelay[i][j]

        elif i==113:# S3 (East)  # Green fiber @ pod S13 is reference channel.
            posDet[i]=[105.1,    395.7,    2632.4]
            podN[i]= [153, 0,0]
            machines[i]= [113,0,0]
            cable[i]= [+80,0,0 ] #R6874
            fiberdelay[i]= [0,0,0] #Reference   
            isSci[i]= 1
            delay[i]=[0,0,0]
            j=int(0)
            for j, elt in enumerate(delay[i]):
                delay[i][j]= cable[i][j]+fiberdelay[i][j]

        elif i==114:# S4 (West)
            posDet[i]=[-74.4,    585.3,    2639.1]
            podN[i]= [147,0,0]
            machines[i]= [114,0,0]
            cable[i]= [+89,0,0] #R6943
            fiberdelay[i]= [-298,0,0] # R6920
            isSci[i]= 1
            delay[i]=[0,0,0]
            j=int(0)
            for j, elt in enumerate(delay[i]):
                delay[i][j]= cable[i][j]+fiberdelay[i][j]

        elif i==115:# S5 (Center)
            posDet[i]=[-4.2,    595.3,    2638.5]
            podN[i]= [147,0,0]
            machines[i]= [115,0,0]
            cable[i]= [+62,0,0 ] #R6944
            fiberdelay[i]= [-297,0,0] # R6923
            isSci[i]= 1
            delay[i]=[0,0,0]
            j=int(0)
            for j, elt in enumerate(delay[i]):
                delay[i][j]= cable[i][j]+fiberdelay[i][j]

        elif i==116:# S6 (East)
            posDet[i]=[100.1,    581.8,    2637.1]
            podN[i]= [147, 0,0]
            machines[i]= [116,0,0]
            cable[i]= [+112,0,0 ]  #R6945
            fiberdelay[i]= [-296,0,0] # R6919   
            isSci[i]= 1
            delay[i]=[0,0,0]
            j=int(0)
            for j, elt in enumerate(delay[i]):
                delay[i][j]= cable[i][j]+fiberdelay[i][j]

        else:
            print('SetupCharacteristics error : Unknow detector')


    return posDet, podN, machines, delay, cable, isSci


if __name__ == "__main__": #si le module nest pas importe mais execute seul
    print(RunSetupBuilder(7005))

