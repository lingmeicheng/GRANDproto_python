#!/home/sandra/anaconda3/bin/ipython

from SharedGlobals import *
from ToolsFunctions import *
import pickle
import numpy as np
import os.path
import struct #Interpret strings as packed binary data
from math import *
from scipy import stats
import decimal
import matplotlib.pyplot as plt


def CoinctableFileBuilder(nrun): #mettre Struct en argument au lieu de nrun quand termine

    #pour ouvrir struct au lieu de passer par DstBuilderLauncher.py

    runstr='R00'+str(nrun)
    dstfile="{0}dst_{1}".format(DST_DATA_PATH,runstr)
    with open(dstfile,'rb') as fs:
        unpick=pickle.Unpickler(fs)
        Struct=unpick.load()

    #load of dictionnaries

    det=Struct['Setup']['Det'] #det is a list of ndets dictionnaries

    indant=[]
    indscint=[]
    for i,elt in enumerate(det):
        if elt['isScint']==0:
            indant.append(i)
        else:
            indscint.append(i)
    indant=np.asarray(indant)
    indscint=np.asarray(indscint)
    print(indant)

    IdCoinc=Struct['Coinc']['IdCoinc']
    tag=Struct['Coinc']['Det']['Tag']
    uTime=Struct['Coinc']['Det']['UnixTime']
    DetId=Struct['Coinc']['Det']['Id']
    TrigTime=Struct['Coinc']['Det']['TrigTime']
    EvtId=Struct['Coinc']['Det']['Evt']
    #CoefCor=Struct['Coinc']['Det']['CoefCor'] #need to be implemented
    #TrigCor=Struct['Coinc']['Det']['TrigCor'] #need to be implemented
    Sat=Struct['Coinc']['Det']['Sat']
    AmpMax=Struct['Coinc']['Det']['AmpMax']


    #create file with antenna coordinates

    if AnalysisType==0:
        txtfile="{0}coord_antennas.txt".format(TEXT_PATH)
    elif AnalysisType==1:
        txtfile="{0}coord_detectors.txt".format(TEXT_PATH)
    elif AnalysisType==2:
        txtfile="{0}coord_scints.txt".format(TEXT_PATH)
    
    with open (txtfile, 'w') as file :
        for i,elt in enumerate(det):
            file.write(str(elt['Id'])+' '+str(elt['X'])+' '+str(elt['Y'])+' '+str(elt['Z'])+'\n')


    #create coinctable file
    thresh=4
    runstr2='R'+str(nrun)
    txtfile="{0}{1}_coinctable.txt".format(TEXT_PATH,runstr2)
    with open (txtfile, 'w') as file :
        for i in range(0,np.shape(IdCoinc)[0]):

            if AnalysisType==0:
                iDetsIn=indant[tag[i,indant]==1]
            elif AnalysisType==1:
                iDetsIn=np.nonzero(tag[i,:]==1)[0]
            elif AnalysisType==2:
                iDetsIn=indscint[tag[i,indscint]==1]
            #print(iDetsIn)

            if np.shape(iDetsIn)[0]>=thresh:
                indtag=np.nonzero(tag[i,:]==1)[0]
                #print('indtag',indtag)

                correlTime=0
                correlCoef=0

                for j in range(0,np.shape(indtag)[0]):
                    chevt=np.nonzero(EvtId[i,indtag[j],:])[0]
                    #print(chevt,EvtId[i,indtag[j],chevt],TrigTime[i,indtag[j],chevt])
                    maxtrig=max(TrigTime[i,indtag[j],chevt])
                    indch=np.nonzero(TrigTime[i,indtag[j],chevt]==maxtrig)[0][0]
                    indtrig=chevt[indch]
                    #print(maxtrig,indtrig,TrigTime[i,indtag[j],indtrig])
                    #maxtrig=maxtrig-mintrigevent
                    #print(i,indtag[j],indtrig,np.shape(indtag)[0])
                    if CORREL:
                        correlTime=TrigCor[i,indtag[j],indtrig]
                        correlCoef=CoefCor[i,indtag[j],indtrig]


                    trigtimedec = decimal.Decimal(str(TrigTime[i,indtag[j],indtrig]))
                    ampmaxdec = decimal.Decimal(str(AmpMax[i,indtag[j],indtrig]))
                         
                    trigtimeround = trigtimedec.quantize(decimal.Decimal('.01'), rounding=decimal.ROUND_HALF_UP)
                    ampmaxround = ampmaxdec.quantize(decimal.Decimal('.01'), rounding=decimal.ROUND_HALF_UP)


                    file.write(str(uTime[i,indtag[j]])+' '+str(DetId[i,indtag[j]])+' '+str(EvtId[i,indtag[j],indtrig])+' '+str(IdCoinc[i])+' '+str(trigtimeround)+' '+str(correlTime)+' '+str(correlCoef)+' '+str(Sat[i,indtag[j],indtrig])+' '+str(ampmaxround)+ ' 0 0 0\n')





    return Struct



def ReconsLoader(nrun): #mettre Struct en argument au lieu de nrun quand termine

    #pour ouvrir struct au lieu de passer par DstBuilderLauncher.py
    runstr='R00'+str(nrun)
    processfile="{0}dst_{1}".format(DST_DATA_PATH,runstr)
    with open(processfile,'rb') as fs:
        unpick=pickle.Unpickler(fs)
        Struct=unpick.load()


    runstr2='R'+str(nrun)
    errort=ErrorTrig/FSAMPLING*1e9

    #load of dictionnaries

    IdCoinc=Struct['Coinc']['IdCoinc']
    det=Struct['Setup']['Det'] #det is a list of ndets dictionnaries
    EvtId=Struct['Coinc']['Det']['Evt']    
    TrigTime=Struct['Coinc']['Det']['TrigTime'] #need to be implemented
    X=[]
    Y=[]
    Z=[]

    for j in range(0,3):
        for i,elt in enumerate(det):
            X.append(elt['X'])
            Y.append(elt['Y'])
            Z.append(elt['Z'])
    X=np.asarray(X)
    Y=np.asarray(Y)
    Z=np.asarray(Z)

    #initialisation of arrays

    Flag=np.zeros(np.shape(IdCoinc)[0])
    Timep=np.zeros(np.shape(IdCoinc)[0])-1
    Multp=np.zeros(np.shape(IdCoinc)[0])-1
    Thetap=np.zeros(np.shape(IdCoinc)[0])-1
    Dthetap=np.zeros(np.shape(IdCoinc)[0])-1
    Phip=np.zeros(np.shape(IdCoinc)[0])-1
    Dphip=np.zeros(np.shape(IdCoinc)[0])-1
    Chi2=np.zeros(np.shape(IdCoinc)[0])-1
    Signif=np.zeros(np.shape(IdCoinc)[0])-1

    chi2ndfp=np.zeros(np.shape(IdCoinc)[0])-1
    slopep=np.zeros(np.shape(IdCoinc)[0])-1

    if AnalysisType==0:
        Flags=np.zeros(np.shape(IdCoinc)[0])
        Times=np.zeros(np.shape(IdCoinc)[0])-1
        Mults=np.zeros(np.shape(IdCoinc)[0])-1
        X0=np.zeros(np.shape(IdCoinc)[0])-1
        Y0=np.zeros(np.shape(IdCoinc)[0])-1
        Z0=np.zeros(np.shape(IdCoinc)[0])-1
        Rhos=np.zeros(np.shape(IdCoinc)[0])-1
        Thetas=np.zeros(np.shape(IdCoinc)[0])-1
        Phis=np.zeros(np.shape(IdCoinc)[0])-1
        Chi2s=np.zeros(np.shape(IdCoinc)[0])-1

        chi2ndfs=np.zeros(np.shape(IdCoinc)[0])-1
        slopes=np.zeros(np.shape(IdCoinc)[0])-1


    #load of plane recons

    if AnalysisType==0:
        txtfile="{0}{1}_planerecons.txt".format(TEXT_PATH,runstr2)
    if AnalysisType==1:
        txtfile="{0}{1}_planereconsHybrid.txt".format(TEXT_PATH,runstr2)
    if AnalysisType==2:
        txtfile="{0}{1}_planereconsScint.txt".format(TEXT_PATH,runstr2)

    with open (txtfile, 'r') as file :
        txt = file.read()
    txt=txt.split('\n')
    
    idcoinc=[]   
    timep=[]
    multp=[]
    thetap=[]
    dthetap=[]
    phip=[]
    dphip=[]
    chi2=[]
    signif=[]
    print(txt[0])
    for l in range(0,len(txt)-1):
        line=txt[l].split(' ')
        linegoods=[]
        for k,elt in enumerate(line):
            if elt!='':
                linegoods.append(elt)
        #print(linegoods)
        idcoinc.append(int(linegoods[0]))
        timep.append(float(linegoods[1]))
        multp.append(int(linegoods[2]))
        thetap.append(float(linegoods[3]))
        dthetap.append(float(linegoods[4]))
        phip.append(float(linegoods[5]))
        dphip.append(float(linegoods[6]))
        chi2.append(float(linegoods[7]))
        signif.append(float(linegoods[8]))



    for i in range(0,len(idcoinc)):
        indcoinc=np.nonzero(IdCoinc==idcoinc[i])[0][0]
        print(indcoinc)
        Flag[indcoinc]=1
        Timep[indcoinc]=timep[i]
        Multp[indcoinc]=multp[i]
        Thetap[indcoinc]=thetap[i]
        Dthetap[indcoinc]=dthetap[i]
        Phip[indcoinc]=phip[i]
        Dphip[indcoinc]=dphip[i]
        Chi2[indcoinc]=chi2[i]
        Signif[indcoinc]=signif[i]


        #transform the direction of the wave in the direction to the source and cancel degeneresence
        if Thetap[indcoinc]>90:
            Thetap[indcoinc]=180-Thetap[indcoinc]

        if Phip[indcoinc]>=180:
            Phip[indcoinc]=Phip[indcoinc]-180
        else:
            Phip[indcoinc]=Phip[indcoinc]+180

        evtflat=EvtId[indcoinc,:,:].flatten(1)
        Xcoinc=X[evtflat>0]
        Ycoinc=Y[evtflat>0]
        Zcoinc=Z[evtflat>0]
        detPos=np.zeros((len(Xcoinc),3))
        detPos[:,0]=Xcoinc
        detPos[:,1]=Ycoinc
        detPos[:,2]=Zcoinc
        print(detPos)
        k=np.asarray([sin(Phip[indcoinc]*np.pi/180)*sin(Thetap[indcoinc]*np.pi/180),  -cos(Phip[indcoinc]*np.pi/180)*sin(Thetap[indcoinc]*np.pi/180),  -cos(Thetap[indcoinc]*np.pi/180)])
        print(k)
        plandelays=np.dot(detPos,k)

        plandelays=plandelays-min(plandelays)
        plandelays=plandelays/C0*1e9
        print(plandelays)
        if CORREL:
            trigcorflat=TrigCor[indcoinc,:,:].flatten(1)
            delays=trigcorflat[evtflat>0]
        else:
            trigtimeflat=TrigTime[indcoinc,:,:].flatten(1)  
            delays=trigtimeflat[evtflat>0]            
        delays=delays/FSAMPLING*1e9; 
        delays=delays-min(delays);
        print(delays)
        #fit=stats.linregress(plandelays,delays)
        fit=np.linalg.lstsq(plandelays[:,np.newaxis],delays)   
        #slopep[indcoinc]=fit[0]
        slopep[indcoinc]=fit[0][0]
        print('slope=',slopep[indcoinc])
        if 0:
            x=np.linspace(min(plandelays),max(plandelays))
            plt.plot(plandelays,delays,'*',x,x*slopep[indcoinc])
            #plt.plot(x,x*slopep)
            plt.show()
        chi2ndfp[indcoinc]=sum((delays-plandelays)**2)/(errort**2)/(np.shape(delays)[0]-1);


    PlanRecons=dict()

    PlanRecons['Flag']=Flag
    PlanRecons['L']=Multp
    PlanRecons['T']=Timep
    PlanRecons['Theta']=Thetap
    PlanRecons['dTheta']=Thetap
    PlanRecons['Phi']=Phip
    PlanRecons['dPhi']=Dphip
    PlanRecons['Chi2']=Chi2
    PlanRecons['Signif']=Signif
    PlanRecons['Chi2Delay']=chi2ndfp
    PlanRecons['SlopeDelay']=slopep

    #load of spherical recons for antennas only

    if AnalysisType==0:
        txtfile="{0}{1}_sphrecons.txt".format(TEXT_PATH,runstr2)
        with open (txtfile, 'r') as file :
            txt = file.read()
        txt=txt.split('\n')

        idcoinc=[]
        times=[]
        mults=[]
        x0=[]
        y0=[]
        z0=[]
        chi2s=[]

        for l in range(0,len(txt)-1):
            line=txt[l].split(' ')
            linegoods=[]
            for k,elt in enumerate(line):
                if elt!='':
                    linegoods.append(elt)
            #print(linegoods)
            idcoinc.append(int(linegoods[0]))
            times.append(float(linegoods[1]))
            mults.append(int(linegoods[2]))
            x0.append(float(linegoods[3]))
            y0.append(float(linegoods[4]))
            z0.append(float(linegoods[5]))
            chi2s.append(float(linegoods[7]))

        for i in range(0,len(idcoinc)):
            indcoinc=np.nonzero(IdCoinc==idcoinc[i])[0][0]

            Flags[indcoinc]=1
            Times[indcoinc]=times[i]
            Mults[indcoinc]=mults[i]
            X0[indcoinc]=x0[i]
            Y0[indcoinc]=y0[i]
            Z0[indcoinc]=z0[i]
            Chi2s[indcoinc]=chi2s[i]

            if Z0[indcoinc]<REFALT:
                Z0[indcoinc]=REFALT+(REFALT-Z0[indcoinc]) #because of the ground-reflectional symmetry in recons

            Rhos,Thetas,Phis=Convert2Sph(X0[indcoinc],Y0[indcoinc],Z0[indcoinc])
            #Xs=[X0[indcoinc],Y0[indcoinc],Z0[indcoinc]]

            
            evtflat=EvtId[indcoinc,:,:].flatten(1)
            Xcoinc=X[evtflat>0]
            Ycoinc=Y[evtflat>0]
            Zcoinc=Z[evtflat>0]
            detPos=np.zeros((len(Xcoinc),3))
            detPos[:,0]=Xcoinc
            detPos[:,1]=Ycoinc
            detPos[:,2]=Zcoinc

            #oneM=np.zeros(np.shape(Xcoinc)[0])+1
            #XsM=np.dot(oneM,Xs)
            Xs=np.zeros((len(Xcoinc),3))
            Xs[:,0]=X0[indcoinc]
            Xs[:,1]=Y0[indcoinc]
            Xs[:,2]=Z0[indcoinc]
            #print(detPos)            
            sphdelays=np.sqrt(np.sum( (detPos-Xs)**2, axis=1 ))
            sphdelays=sphdelays-min(sphdelays)
            sphdelays=sphdelays/C0*1e9
            #print(sphdelays)

            if CORREL:
                trigcorflat=TrigCor[indcoinc,:,:].flatten(1)
                delays=trigcorflat[evtflat>0]
            else:
                trigtimeflat=TrigTime[indcoinc,:,:].flatten(1)  
                delays=trigtimeflat[evtflat>0]            
            delays=delays/FSAMPLING*1e9; 
            delays=delays-min(delays);

            #fit=stats.linregress(sphdelays,delays)        
            fit=np.linalg.lstsq(sphdelays[:,np.newaxis],delays) 
            #slopes[indcoinc]=fit[0]
            slopes[indcoinc]=fit[0][0]
            chi2ndfs[indcoinc]=sum((delays-sphdelays)**2)/(errort**2)/(np.shape(delays)[0]-1);


        SphRecons=dict()

        SphRecons['Flag']=Flags
        SphRecons['L']=Mults
        SphRecons['T']=Times
        SphRecons['Rho']=Rhos
        SphRecons['Theta']=Thetas
        SphRecons['Phi']=Phis
        SphRecons['X0']=X0
        SphRecons['Y0']=Y0
        SphRecons['Z0']=Z0
        SphRecons['Chi2']=Chi2
        SphRecons['Chi2Delay']=chi2ndfs
        SphRecons['SlopeDelay']=slopes


    #etablishment of dictionnaries hierarchy

    if AnalysisType==0:
        Radio=dict()
        Radio['Radio']=PlanRecons
        Struct['Coinc']['PlanRecons']=Radio
        Struct['Coinc']['SphRecons']=SphRecons

    elif AnalysisType==1:
        Hybrid=dict()
        Hybrid['Hybrid']=PlanRecons
        Struct['Coinc']['PlanRecons']=Hybrid

    elif AnalysisType==2:
        Struct['Coinc']['PlanRecons']=PlanRecons


    #saving the struct dictionnary in dstfile

    dstfile="{0}dst_{1}".format(DST_DATA_PATH,runstr)
    with open(dstfile,'wb') as fs:
        pick=pickle.Pickler(fs)
        pick.dump(Struct)


    return Struct


if __name__ == "__main__": #si le module nest pas importe mais execute seul
    CoinctableFileBuilder(7005)
    #ReconsLoader(7005)



