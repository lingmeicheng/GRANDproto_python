# Script to analyze GP35 slow control data 
# OMH Aug. 29, 2018
# Updated Nov 16, on same structure as minBias analysis. Now analysing one board at a time only + reduce info to one txt file.

import os
from os.path import expanduser
import time
import datetime
#import mx.DateTime
import sys
import math
import numpy as np
import pylab as pl
import matplotlib.style
import matplotlib as mpl
mpl.style.use('classic')

def loopSLCEvents(boardID,RUNID):
   #datadir = "/home/pastsoft/data/"
   datadir = "/home/martineau/GRAND/GRANDproto35/data/ulastai/"
   filename = datadir+"S"+RUNID+"_b"+boardID+".data.txt"   # To be modified
   if os.path.isfile(filename) is False:
     print 'File ',filename,'does not exist. Aborting.'
     return

   # Load data file
   print 'Scanning slow control datafile',filename
   with open(filename,"r") as f:
   	   evts = f.read().split('-----------------')
   nevts=len(evts)-1
   
   # Initializing arrays
   date = []
   #np.zeros(shape=(np.size(evts)),dtype = np.int32)
   utcsec = []
   VPower = []
   Temp = []
   maxCoarse = []
   TrigRate = []
   Th = np.zeros(shape=(nevts,3,2))
   
   ## Loop on all events; filing in arrays
   #j = 0
   for i in range(nevts):
      evt = evts[i]
      if len(evt)>60:  # Skip ACK words (size=59)
           evtsplit = evt.split('\n')
	   IP = evtsplit[2][3:]
	   board = int(IP[-2:]);
	   if board != int(boardID):
	     #print 'This is board {0}, skiping it (analysing board {1} only)'.format(board,boardID)
	     continue
		   
	   # Now reducing data
	   # Now build time info... Alternative: mx.DateTime.DateTimeFrom(date[0])
	   date=evtsplit[1]
           thisDatetime = datetime.datetime.strptime(date, '%a %b %d %H:%M:%S %Y GMT')  # Build DateTime object
           utcsec.append(time.mktime(thisDatetime.timetuple()))  # Build UTC second
	   
	   power = np.zeros(shape=(6,))
	   for k in range(0,6):
             power[k]=evtsplit[k+3].split(':')[1]
           VPower.append(power)
	   
	   #for k in range(0,3):
           #  thsplit = evtsplit[k+9].split(':')
           #  Th[j,k,0] = thsplit[1].split(' ')[0]
           #  Th[j,k,1] = thsplit[1].split(' ')[1]

           Temp.append(evtsplit[12].split(':')[1])

           trate = np.zeros(shape=(7,))
	   for k in range(0,7):  # 7 columns: overall rate + 6 individual channels
             trate[k] = evtsplit[13+k].split(':')[1]
           TrigRate.append(trate)
	   
	   maxCoarse.append(evtsplit[20].split(':')[1])

           #j = j+1
	   
	   
   #nevts = j
   #date = np.array(date)
   # Build proper temperature info from raw temperature data... Already done in DAQ!!!
   #Temp = np.asarray(Temp)
   #hraw = [hex(int(a)) for a in Temp]  # Transfer back to hexadecimal
   #braw = [bin(int(a)) for a in Temp]
   #draw = [twos_comp(int(a,16), 13) for a in hraw] #2s complements
   
   # Write to file
   utcsec = np.asarray(utcsec,dtype=int)
   TrigRate = np.asarray(TrigRate,dtype=float)
   Temp = np.asarray(Temp,dtype=float)
   VPower = np.asarray(VPower,dtype=float)
   maxCoarse = np.asarray(maxCoarse,dtype=int)
   nev = len(utcsec)
   conc = np.concatenate((utcsec.reshape(nev,1),Temp.reshape(nev,1),VPower.reshape(nev,6),TrigRate.reshape(nev,7),maxCoarse.reshape(nev,1),),axis=1)   # Concatenate results
   #conc = conc.reshape(np.size(utcsec),16) # 
   filename = 'SLC_b{0}.txt'.format(boardID)
   np.savetxt(filename,conc,fmt='%3.2f')  # Write to file

def displaySLC(boardID):
   home = expanduser("~")
   resdir = home+"/GRAND/GRANDproto35/data/ulastai/"
   resfile = resdir+"SLC_b"+str(boardID)+".txt"
   print "Calling displaySLC(). Will display SLC result file {0}".format(resfile)
   
   # Load data
   #open(resfile,'ab')
   a = np.loadtxt(resfile)
   utc = a[:,0]
   temp = a[:,1]
   V = a[:,2:8]
   trig = a[:,8:16]
   mCoarse = a[:,15]
   
   print temp
   print V
   print trig[:,0]
   print mCoarse
         
   #time = (time-min(time[time>0]))/60
   Triglabel = ['Total','Ch1+','Ch2+','Ch3+','Ch1-','Ch2-','Ch3-']
   Voltlabel = ['Main','-3V','+4V','LNA1','LNA2','LNA3']
   print 'SLC info for board',boardID,' in period:'
   print datetime.datetime.fromtimestamp(min(utc)).strftime('%y/%m/%d - %H:%M:%S UTC'),' to ',datetime.datetime.fromtimestamp(max(utc)).strftime('%y/%m/%d - %H:%M:%S UTC')

   # Time ticks
   nticks = 8
   ind = np.linspace(min(utc),max(utc),nticks)
   date = [datetime.datetime.fromtimestamp(ux).strftime('%H:%M') for ux in ind]

   # Temperature plot
   pl.figure(1)
   pl.plot(utc,temp)
   pl.xticks(np.linspace(min(utc),max(utc),nticks), date)
   pl.xlim(min(utc)-1,max(utc)+1)
   pl.grid(True)
   pl.legend(loc='best')
   pl.xlabel('Time')
   pl.ylabel('Temperature ($^{\circ}$C)')
   pl.title('Board temperature')
   pl.savefig('temp.png')	

   # Voltage plot
   pl.figure(2)
   for i in range(np.shape(V)[1]):
        sub=321+i
   	pl.subplot(sub)
	pl.plot(utc,V[:,i],lw=2,label=Voltlabel[i])
	pl.grid(True)
	if i>3:
  	  pl.xlabel('Time')
	pl.xticks(np.linspace(min(utc),max(utc),nticks), date)
	pl.title(Voltlabel[i])
        pl.ylabel('Voltage (V)')
   pl.savefig('voltage.png')	

   pl.figure(3)  #Trig Rate
   # Plotting total trig rate only. 
   # pl.plot(time[sel],TrigRate[sel,0][0],lw=2,label=id) 
   # Plotting all individual trig rates. 
   labs = ['Total','X+','Y+','Z+','X-','Y-','Z-']
   for ch in range(7):
     pl.plot(utc,trig[:,ch],lw=2,label=labs[ch])

   pl.xticks(np.linspace(min(utc),max(utc),nticks), date)
   pl.grid(True)
   pl.ylabel('Total trig rate (Hz)')
   pl.xlabel('Time')
   pl.legend(loc='best')
   pl.savefig('trig.png')	

   pl.draw()
   pl.show()

   return

def twos_comp(val, bits):
   """compute the 2's compliment of int value val"""
   if (val & (1 << (bits - 1))) != 0: # if sign bit is set e.g., 8bit: 128-255
     val = val - (1 << bits)	    # compute negative value
   return val


if __name__ == '__main__':
#     if len(sys.argv)!=3:
#       print "Usage: >loopSLCEvents RUNID BOARDID"
#     else:  
       #loopSLCEvents(sys.argv[1],sys.argv[2])
       displaySLC(sys.argv[1])
