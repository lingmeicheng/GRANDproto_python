# Script to analyze GP35 slow control data 
# OMH Aug. 29, 2018

import os
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

def loop(RUNID,boardID):
   datadir = "/home/pastsoft/data/"
   #datadir = "/home/martineau/GRAND/GRANDproto35/data/ulastai/"
   filename = datadir+"S"+RUNID+"_b06.data.txt"   # To be modified
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
   board = np.zeros(shape=(np.size(evts)),dtype = np.int32)
   time = np.zeros(shape=(np.size(evts)))
   VPower = np.zeros(shape=(nevts,6))
   Th = np.zeros(shape=(nevts,3,2))
   Temp = np.zeros(shape=(np.size(evts)))
   maxCoarse = np.zeros(shape=(np.size(evts)))
   TrigRate = np.zeros(shape=(nevts,7))
   
   ## Loop on all events; filing in arrays
   j = 0
   for i in range(j,nevts):
      evt = evts[i]
      if len(evt)>60:  # Skip ACK words (size=59)
           evtsplit = evt.split('\n')
	   thisDate=evtsplit[1]
           date.append(thisDate)
	   # Now build time info... Alternative: mx.DateTime.DateTimeFrom(date[0])
	   thisDatetime = datetime.datetime.strptime(thisDate, '%a %b %d %H:%M:%S %Y GMT')
	   time[j] = (thisDatetime-datetime.datetime(2018,8,27)).total_seconds()
           IP = evtsplit[2].split(':')[1]
           board[j] = int(IP[-2:]);

	   for k in range(0,6):
             VPower[j,k]=evtsplit[k+3].split(':')[1]

	   for k in range(0,3):
             thsplit = evtsplit[k+9].split(':')
             Th[j,k,0] = thsplit[1].split(' ')[0]
             Th[j,k,1] = thsplit[1].split(' ')[1]

           Temp[j] = evtsplit[12].split(':')[1]

           for k in range(0,7):  # 7 columns: overall rate + 6 individual channels
             TrigRate[j,k] = evtsplit[13+k].split(':')[1]

	   maxCoarse[j] = evtsplit[20].split(':')[1]

           j = j+1
	   
   nevts = j
   date = np.array(date)
   # Build proper temperature info from raw temperature data
   hraw = [hex(int(a)) for a in Temp]  # Transfer back to hexadecimal
   braw = [bin(int(a)) for a in Temp]
   draw = [twos_comp(int(a,16), 13) for a in hraw] #2s complements

   time = (time-min(time[time>0]))/60
   Triglabel = ['Total','Ch1+','Ch2+','Ch3+','Ch1-','Ch2-','Ch3-']
   Voltlabel = ['Main','-3V','+4V','LNA1','LNA2','LNA3']
   boards = set(board[np.where(board>0)])
   print 'Boards in run:',list(boards)
   print 'Run start:', date[0]
   
   # Loop on all boards
   #for id in boards:
   for id in [int(boardID)]:
   	   sel = np.where(board == id)
	   date_end = date[sel[0][-1]]
   	   print 'Run stop:',date_end,'for board',id,' (',np.size(sel),'measurements)'

	   # Temperature plot
	   pl.figure(1)
   	   pl.plot(time[sel],Temp[sel],label=id)
	   pl.grid(True)
	   pl.legend(loc='best')
	   pl.xlabel('Time (mn)')
	   pl.ylabel('Temperature ($^{\circ}$C)')
	   pl.title('Board temperature')
	   pl.savefig('temp.png')	

	   # Voltage plot
	   pl.figure(2)
	   for i in range(np.shape(VPower)[1]):
	        sub=321+i
	   	pl.subplot(sub)
		pl.plot(time[sel],VPower[sel,i][0],lw=2,label=id)
		pl.grid(True)
		if i ==0:
		  pl.legend(loc='best')
		if i>3:
  		  pl.xlabel('Time (mn)')
		pl.title(Voltlabel[i])
 	        pl.ylabel('Voltage (V)')
	   pl.savefig('voltage.png')	

	   pl.figure(3)  #Trig Rate
           # Plotting total trig rate only. 
	   # pl.plot(time[sel],TrigRate[sel,0][0],lw=2,label=id) 
	   # Plotting all individual trig rates. 
           labs = ['Total','X+','Y+','Z+','X-','Y-','Z-']
	   for ch in range(7):
	     pl.plot(time[sel],TrigRate[sel,ch][0],lw=2,label=labs[ch])
	   
	   #pl.yscale('log')
	   pl.grid(True)
	   pl.ylabel('Total trig rate (Hz)')
	   pl.xlabel('Time (mn)')
	   pl.legend(loc='best')
	   #pl.title('Board '+str(id))
	   pl.savefig('trig.png')	
	   
	   thisMaxCoarse = maxCoarse[sel]
	   # MaxCoarse info
	   pl.figure(4)  #MaxCoarse
	   pl.subplot(211)
	   pl.plot(time[thisMaxCoarse>0],thisMaxCoarse[thisMaxCoarse>0])
	   pl.ylabel('MaxCoarse')
	   pl.xlabel('Time (mn)')
	   pl.subplot(212)
	   pl.hist(thisMaxCoarse[thisMaxCoarse>0])
	   pl.xlabel('MaxCoarse')
	   #print 'maxCoarse counter = ',np.mean(thisMaxCoarse[thisMaxCoarse>0]),'+-',np.std(thisMaxCoarse[thisMaxCoarse>0]),' vs 1.25e+8'

   pl.draw()
   pl.show()

   return

def twos_comp(val, bits):
   """compute the 2's compliment of int value val"""
   if (val & (1 << (bits - 1))) != 0: # if sign bit is set e.g., 8bit: 128-255
     val = val - (1 << bits)	    # compute negative value
   return val


if __name__ == '__main__':
     if len(sys.argv)!=3:
       print "Usage: >loop RUNID BOARDID"
     else:  
       loop(sys.argv[1],sys.argv[2])
