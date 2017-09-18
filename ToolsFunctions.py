#!/home/sandra/anaconda3/bin/ipython

from math import *
from SharedGlobals import *
from scipy.fftpack import rfft, irfft, rfftfreq
import numpy as np
import scipy.signal

def UnixSecs2Date(unixtime):

    return 0

#############

def Convert2Sph(x,y,z):
#convention phi=0 on y, phi=90 on -x
    x=x-REFWE
    y=y-REFSN
    z=z-REFALT

    r=sqrt(x**2+y**2+z**2)
    theta=acos(z/r)
    sinphi=-x/(r*sin(theta))
    if sinphi<0:
        phi=asin(sinphi)+2*pi
    else:
        phi=asin(sinphi)

    theta=theta*180/pi
    phi=phi*180/pi

    return r,theta,phi

#############

def PassBand(vs,ts,fmin,fmax):
    tstep=np.mean(np.diff(ts))
    F=rfftfreq(len(vs))/tstep #len(vs) points between 0 and 1/2tstep=0.5e10Hz
    VS=rfft(vs)
    VS[F<fmin]=0
    VS[F>fmax]=0
    vs=irfft(VS)

    return(vs)

#############

def Agglomerate(x,g):

    g1=g
    k=g+1
    k1=g1+1
    
    #filter isolated 1
    k2=2*k1+1
    if g1>1:
        xx=np.zeros(len(x)+k2)
        xx[0:len(x)]=x
        inter=np.zeros(k2)+1
        xx=scipy.signal.lfilter(inter,1,xx)
        xx=xx[k1:len(xx)-k1-1]
        boolx=np.zeros(len(x),bool)
        for i in range(0,len(x)):
            boolx[i]=(x[i]==1 and xx[i]>1)
        x[:]=0
        x[boolx]=1

    #filter holes
    mask=np.arange(k-1,-1,-1)
    mask=2**mask

    yb=scipy.signal.lfilter(mask,1,x) #left
    K=yb>0
    J=yb==0
    yb[K]=k-1-np.floor(np.log2(yb[K]))
    yb[J]=k+1

    ya=scipy.signal.lfilter(mask,1,x[::-1]) #right
    ya=ya[::-1]
    K=ya>0
    J=ya==0
    ya[K]=k-1-np.floor(np.log2(ya[K]))
    ya[J]=k+1

    y=np.zeros(len(x))
    booly=np.zeros(len(y),bool)
    for i in range(0,len(y)):
        booly[i]= (x[i]==1 or (ya[i]+yb[i])<=k)
    y[booly]=1

    return(y)

#############



if __name__ == "__main__": #si le module nest pas importe mais execute seul
    #test
    x=[1,1,1,1,0,0,0,0,0,1,1,1,0,1,1,1,0,0,0,0,1,1,0,0]
    x=np.asarray(x,float)
    print(x)
    g=3
    y=Agglomerate(x,g) #we can input here as many arguments as needed

    print(y)













