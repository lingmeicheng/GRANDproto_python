#!/usr/bin/env python
import os
import glob
import re

import numpy as np
import pylab as pl


def fetchPath(): 
  # path to the directory (relative or absolute)
  dirpath = "/home/martineau/GRAND/GRANDproto35/data/tests_STAE_final"
  dirList = [];
  print "Walking through directories..."
  for dirname, dirnames, _ in os.walk(dirpath):  # Walk through the structure
      if len(dirnames)>1:
  	continue  # Skip first part where we list all boardid dirs, and jump one level down...
      for subdirname in dirnames:
  	  #print "dir=",dirname
  	  #print "sudir=",subdirname
  	  dirList.append( os.path.join(dirname, subdirname))
	  
      #if len(dirList)>5:
      # break

  return dirList
  
     
def anaThresh(dirList,option):
     print option
     board = np.empty((0,1),int)
     effX = np.empty((0,1),float)
     
     for d in dirList:
       if option == 's':
         threshfile =  glob.glob(d+'/threshs_b*.txt')
       if option == 'n':
         threshfile =  glob.glob(d+'/threshn_b*.txt')
       if len(threshfile)>0:
         threshfile = threshfile[0]
	 print "Reading threshold file",threshfile
	  
         bid = threshfile.split('_b')[1]
         bid = bid.split('.')[0]
         board = np.append(board,int(bid))
	 # Read data
 	 a = np.loadtxt(threshfile)
 	 pat = a[:,0]
 	 th = a[:,1]
 	 nb = a[:,2]
 	 y = np.asarray(nb)/(10*10*2)  #10Hz square
 	 pati = [1000,10000,100000]
         ch = ['X+','Y+','Z+']
	 for i in range(np.size(pati)):
           sel = np.where(pat==pati[i])
	   print ch[i],th[sel],y[sel]
	   if option == 's':
	     effX = np.append(effX,y[sel][0:3]) 
	   if option == 'n':
	     effX = np.append(effX,y[sel][1:4]*20)   
     
     print 'Boards with data:', len(board)
     if option == 's':
       print 'Ratio in [0.95 - 1.05]:',np.sum( (effX>0.95) & (effX<1.05)),'/',len(effX),'=',float(np.sum( (effX>0.95) & (effX<1.05)))/len(effX)
       print 'Ratio<0.95:',np.sum(effX<0.95),'/',len(effX),'=',float(np.sum(effX<0.95))/len(effX)
       print 'Ratio>1.05:',np.sum(effX>1.05),'/',len(effX),'=',float(np.sum(effX>1.05))/len(effX)
     if option == 'n':
       print 'Rate>1:',np.sum(effX>1),'/',len(effX),'=',float(np.sum(effX>1))/len(effX)

     pl.figure(1)
     if option == 's':
        pl.hist(effX[effX<=2],100)
	pl.xlabel('Signal trigger ratio')
     if option == 'n':
        pl.hist(effX[effX<5],200)
	pl.xlabel('Noise trigger rate (Hz)')
     pl.show()
    
       	  
def anaTrig(anaDir):
     board = np.empty((0,1),int)
     freqR = np.empty((0,1),float)
     
     for d in dirList:
       hrFile =  glob.glob(d+'/R35_b*.data')	# This is the correct file ID for this high trig rate file
       if len(hrFile)>0:
         hrFile = hrFile[0]
         print  hrFile,os.stat(hrFile).st_size
         bid = hrFile.split('_b')[1]
	 bid = bid.split('.')[0]
	 board = np.append(board,int(bid))
	   
	 # Read data
         with open(hrFile,"r") as f:
   	   
	   evts = f.read().split('-----------------')
           nevts = len(evts)-1
           print 'Number of events:',nevts
           evt0 = evts[1]
	   evtsplit = evt0.split('\n')
	   s0=int(evtsplit[6][4:])  # Elapsed seconds since start
	   ts2=int(evtsplit[3][4:])
	   t0 = s0+ts2*8e-9
	   evtEnd = evts[-2]
	   evtsplit = evtEnd.split('\n')
	   send=int(evtsplit[6][4:])  # Elapsed seconds since start
	   ts2=int(evtsplit[3][4:])
	   tend = send+ts2*8e-9
	   tdaq = tend-t0
	   if tend>0:
	     tdaq = tend-t0
	   else:  # No GPS data
	     tdaq = 5.2
	   freqR = np.append(freqR,nevts/tdaq)
	   
     bad = freqR<2500	   
     print 'Boards with evt rate<2.5kHz:',board[bad]
     print len(freqR[freqR>=2500])
     pl.figure(1)
     pl.hist(freqR[freqR>=2500],20)
     pl.xlabel('Event rate (Hz)')
     pl.show()
	   
def anaFreq(dirList):     
     board = np.empty((0,1),int)
     mx = np.empty((0,1),float)
     nonflatall = np.empty((0,3),float)
     dropinfall = np.empty((0,3),float)
     dropsupall = np.empty((0,3),float)
     
     freq = np.arange(20,115,5)  # taken from testfreq.sh
     flat = np.where( (freq>35) & (freq<95) )
     
     for d in dirList:
        freqfile = glob.glob(d+'/freq_b*.txt')	
	if  len(freqfile)>0: # summary file was found 
	  freqfile = freqfile[0]
	  print "Reading frequency file",freqfile
	  bid = freqfile.split('_b')[1]
	  bid = bid.split('.')[0]
	  board = np.append(board,int(bid))
	  
	  m,std,c = matchPatternCalib(freqfile,' = ')
	  m = m.reshape(len(m)/3,3)
	  
	  nonflat = np.empty((1,3),float)
	  dropinf = np.empty((1,3),float)
	  dropsup = np.empty((1,3),float)
	  for i in range(3):
	    mi = m[:,i]
	    nonflat[0,i] = np.std(mi[flat])/np.mean(mi[flat])  # Check divergance to perfect flatness
	    dropinf[0,i] = np.mean(mi[flat])-mi[freq==20]  # Check drop on LF side
	    dropsup[0,i] = np.mean(mi[flat])-mi[freq==110]  # Check drop on HF side
	  nonflatall = np.append(nonflatall,nonflat)  
	  dropinfall = np.append(dropinfall,dropinf)
          dropsupall = np.append(dropsupall,dropsup)
	    
	else:
          print "No summary file in",d
     
     #nonflatall = nonflatall.reshape(len(nonflatall)/3,3)
     print 'Boards with freq data:',len(board)
     
     bad = np.where(dropinfall<0.8)
     bad = np.asarray(bad)/3
     print 'Boards with dropinf<0.8V:',board[bad]
     bad = np.where(dropsupall<0.75)
     bad = np.asarray(bad)/3
     print 'Boards with dropsup<0.8V:',board[bad]
   
     pl.figure(1)
     pl.hist(nonflatall,50)
     pl.xlabel('Flatness ratio')
     #pl.subplot(311)
     #pl.hist(nonflatall[:,0],50)
     
     pl.figure(2)
     pl.subplot(211)
     pl.hist(dropinfall,50)
     pl.xlabel('$\Delta$V LF (V)')
     pl.subplot(212)
     pl.hist(dropsupall,50)
     pl.xlabel('$\Delta$V HF (V)')

     pl.show()
     


def matchPatternCalib(filename,pattern): 
    cal = np.empty((0,1),float)
    m = np.empty((0,1),float)
    std = np.empty((0,1),float)
    for i, line in enumerate(open(filename)):
          #print "Reading line",line
	  if len(line.split(pattern))==1: # pattern not found
	    continue
	  if len(line.split(pattern))==2:  # Only one pattern
	    c = line.split(pattern)[1]
	    cal = np.append(cal,float(c.split(' ')[0]))
	    #print line
	    #print "cal=",cal
	  else:  # 2 patterns
	    a = line.split(pattern)[1]
	    m = np.append(m,float(a.split(' ')[0]))
	    b = line.split(pattern)[2]
	    std = np.append(std,float(b.split(' ')[0]))
	    #print "m=",m
	    #print "std=",std
	    
    #m = m.reshape(3,len(m)/3)
    #std = std.reshape(3,len(std)/3)
    #print len(m),len(std),len(cal)
    return m,std,cal	  

def anaCalib(dirList):
     
     board = np.empty((0,1),int)
     m30 = np.empty((0,1),float)
     std30 = np.empty((0,1),float)
     m60 = np.empty((0,1),float)
     std60 = np.empty((0,1),float)
     m90 = np.empty((0,1),float)
     std90 = np.empty((0,1),float)
     cal = np.empty((0,1),float)
     
     for d in dirList:
        sumfile = glob.glob(d+'/summary*.txt')	
	if  len(sumfile)>0: # summary file was found 
	  sumfile = sumfile[0]
	  print "Reading summary file",sumfile
	  bid = sumfile.split('_b')[1]
	  bid = bid.split('.')[0]
	  board = np.append(board,int(bid))
	  
	  m,std,c = matchPatternCalib(sumfile,' = ')
	  m30 = np.append(m30,m[0:3])
	  std30 = np.append(std30,std[0:3])
	  m60 = np.append(m60,m[3:6])
	  std60 = np.append(std60,std[3:6])
	  m90 = np.append(m90,m[6:9])
	  std90 = np.append(std90,std[6:9])
	  print c	  
	  cal = np.append(cal,c)
        else:
          print "No summary file in",d
	  
     print 'Boards with SLC data:',len(board)
     print 'cal=',len(cal)
     print 'm30=',len(m30)
     print 'm60=',len(m60)
     print 'm90=',len(m90)
     
     #bad = np.where(m90>0.5)
     #print 'Boards with m90>0.5V:',m90[bad]
     
     pl.figure(1)
     pl.subplot(311)
     pl.hist(m30,50)
     pl.xlabel('<V> for AttIndex = 30+30 (V)')
     pl.subplot(312)
     pl.hist(m60,50)
     pl.xlabel('<V> for AttIndex = 60+60 (V)')
     pl.subplot(313)
     pl.hist(m90,50)
     pl.xlabel('<V> for AttIndex = 90+90 (V)')
     
     pl.figure(2)
     pl.subplot(311)
     pl.hist(std30,50)
     pl.xlabel('$\sigma_{V}$ for AttIndex = 30+30 (V)')
     pl.subplot(312)
     pl.hist(std60,50)
     pl.xlabel('$\sigma_{V}$ for AttIndex = 60+60 (V)')
     pl.subplot(313)
     pl.hist(std90,50)
     pl.xlabel('$\sigma_{V}$ for AttIndex = 90+90 (V)')
     
     pl.figure(3)
     pl.hist(cal)
     pl.xlabel('<V$_{pp}$> for calibrator (V)')

     attInd = np.array([30,60,90])
     attdB = -(-0.25*attInd+31.75+3.3)*2
     print attdB
     vin = 0.262/2*pow(10,attdB/20)  #mV
     
     fig = pl.figure(4)
     ax = fig.add_subplot(2, 1, 1)
     ax.plot(vin[0]*np.ones((len(m30),1)),m30,'+')
     ax.plot(vin[1]*np.ones((len(m30),1)),m60,'+')
     ax.plot(vin[2]*np.ones((len(m30),1)),m90,'+')
     pl.grid(True)
     ax.set_xscale('log')
     ax.set_xlabel('$V_{in}$ (V)')
     ax.set_ylabel('$V_{out}$ (V)')
     
     ax = fig.add_subplot(2, 1, 2)
     ax.plot(vin[0]*np.ones((len(m30),1)),std30,'+')
     ax.plot(vin[1]*np.ones((len(m30),1)),std60,'+')
     ax.plot(vin[2]*np.ones((len(m30),1)),std90,'+')
     pl.grid(True)
     ax.set_xscale('log')
     ax.set_xlabel('$V_{in}$ (V)')
     ax.set_ylabel('$\sigma_{out}$ (V)')
     
     m30 = m30.reshape(len(m30)/3,3)
     cal = cal.reshape(len(cal)/3,3)
     pl.figure(5)
     pl.plot(cal[:,0],m30[:,0],'+')
     pl.plot(cal[:,1],m30[:,1],'+')
     pl.plot(cal[:,2],m30[:,2],'+')
     
     pl.show()
     
def matchPatternSLC(filename,pattern): 
    for i, line in enumerate(open(filename)):
      for match in re.finditer(pattern, line):
        return float(line.split(':'))
	  
def anaSLC(dirList):
     board = np.empty((0,1),int)
     V1 = np.empty((0,1),float)
     V2 = np.empty((0,1),float)
     V3 = np.empty((0,1),float)
     V4 = np.empty((0,1),float)
     V5 = np.empty((0,1),float)
     V6 = np.empty((0,1),float)
     T = np.empty((0,1),float)
     for d in dirList:
        slcfile = glob.glob(d+'/S*')
	
	if  len(slcfile)>0: # Slow control file was found 
	  slcfile = slcfile[0]
	  print "Reading SLC data file",slcfile
 	  bid = slcfile.split('board')[1]
	  bid = bid.split('/')[0]
	  board = np.append(board,int(bid))
	  V1 = np.append(V1,matchPatternSLC(slcfile,'VPower1'))
 	  V2 = np.append(V2,matchPatternSLC(slcfile,'VPower2'))
	  V3 = np.append(V3,matchPatternSLC(slcfile,'VPower3'))
 	  V4 = np.append(V4,matchPatternSLC(slcfile,'VPower4'))
	  V5 = np.append(V5,matchPatternSLC(slcfile,'VPower5'))
 	  V6 = np.append(V6,matchPatternSLC(slcfile,'VPower6'))
	  T = np.append(T,matchPatternSLC(slcfile,'Temp'))
     
     
     print 'Boards with SLC data:',len(board)
     print board
     print V1
     print V2
     print V3
     print V4
     print V5
     print V6
     print T
     
     bad = np.where(V4<10)
     print 'Boards with LNA1<10V:',board[bad],V4[bad]
     
     pl.figure(1)
     pl.subplot(231)
     pl.hist(V1)
     pl.xlabel('Board input DC level (V)')
     pl.subplot(232)
     pl.hist(V2)
     pl.xlabel('Board 3V supply (V)')
     pl.subplot(233)
     pl.hist(V3)
     pl.xlabel('Board -4V supply (V)')
     pl.subplot(234)
     pl.hist(V4)
     pl.xlabel('LNA1 voltage (V)')
     pl.subplot(235)
     pl.hist(V5)
     pl.xlabel('LNA2 voltage (V)')
     pl.subplot(236)
     pl.hist(V6)
     pl.xlabel('LNA3 voltage (V)')
     
     pl.figure(2)
     pl.hist(T)
     pl.xlabel('Board temperature ($^{\circ}$C)')
     
     pl.show()

	    
if __name__ == '__main__':
     #dirList = fetchPath()
     dirList = ['/home/martineau/GRAND/GRANDproto35/data/tests_STAE_final/board01/171207_1332', '/home/martineau/GRAND/GRANDproto35/data/tests_STAE_final/board02/180104_1855', '/home/martineau/GRAND/GRANDproto35/data/tests_STAE_final/board03/171206_1125', '/home/martineau/GRAND/GRANDproto35/data/tests_STAE_final/board04/171206_1622', '/home/martineau/GRAND/GRANDproto35/data/tests_STAE_final/board05/171206_1206', '/home/martineau/GRAND/GRANDproto35/data/tests_STAE_final/board06/171206_1225', '/home/martineau/GRAND/GRANDproto35/data/tests_STAE_final/board07/171206_1639', '/home/martineau/GRAND/GRANDproto35/data/tests_STAE_final/board08/171206_1150', '/home/martineau/GRAND/GRANDproto35/data/tests_STAE_final/board09/171206_1242', '/home/martineau/GRAND/GRANDproto35/data/tests_STAE_final/board10/171206_1330', '/home/martineau/GRAND/GRANDproto35/data/tests_STAE_final/board11/171206_1346', '/home/martineau/GRAND/GRANDproto35/data/tests_STAE_final/board12/171205_1318', '/home/martineau/GRAND/GRANDproto35/data/tests_STAE_final/board13/171206_1300', '/home/martineau/GRAND/GRANDproto35/data/tests_STAE_final/board14/171206_1655', '/home/martineau/GRAND/GRANDproto35/data/tests_STAE_final/board15/171207_0918', '/home/martineau/GRAND/GRANDproto35/data/tests_STAE_final/board16/171206_1403', '/home/martineau/GRAND/GRANDproto35/data/tests_STAE_final/board17/171206_1059', '/home/martineau/GRAND/GRANDproto35/data/tests_STAE_final/board18/171206_1423', '/home/martineau/GRAND/GRANDproto35/data/tests_STAE_final/board19/171206_1444', '/home/martineau/GRAND/GRANDproto35/data/tests_STAE_final/board20/171206_1502', '/home/martineau/GRAND/GRANDproto35/data/tests_STAE_final/board21/171206_1518', '/home/martineau/GRAND/GRANDproto35/data/tests_STAE_final/board22/171207_0937', '/home/martineau/GRAND/GRANDproto35/data/tests_STAE_final/board24/171207_1312', '/home/martineau/GRAND/GRANDproto35/data/tests_STAE_final/board25/171206_1547', '/home/martineau/GRAND/GRANDproto35/data/tests_STAE_final/board26/171207_1226', '/home/martineau/GRAND/GRANDproto35/data/tests_STAE_final/board27/171207_0954', '/home/martineau/GRAND/GRANDproto35/data/tests_STAE_final/board28/171207_1256', '/home/martineau/GRAND/GRANDproto35/data/tests_STAE_final/board29/171207_1114', '/home/martineau/GRAND/GRANDproto35/data/tests_STAE_final/board30/171207_1049', '/home/martineau/GRAND/GRANDproto35/data/tests_STAE_final/board31/171207_1033', '/home/martineau/GRAND/GRANDproto35/data/tests_STAE_final/board32/171207_1016', '/home/martineau/GRAND/GRANDproto35/data/tests_STAE_final/board33/171207_1208', '/home/martineau/GRAND/GRANDproto35/data/tests_STAE_final/board34/171207_1131', '/home/martineau/GRAND/GRANDproto35/data/tests_STAE_final/board35/171207_1152']
     print dirList
     
     #anaSLC(dirList)
     #anaCalib(dirList)  
     #anaFreq(dirList)  
     #anaTrig(dirList)
     #anaThresh(dirList,'s')
     anaThresh(dirList,'n')
