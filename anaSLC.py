import os
import time
import datetime
import mx.DateTime
import sys
import math
import numpy as np
import pylab as pl

datadir = "/home/martineau/GRAND/GRANDproto35/data/tests/"

def loop(filename):
   datafile = datadir+filename
   print 'Scanning',datafile

   with open(datafile,"r") as f:
   	   evts = f.read().split('-----------------')

   nevts=len(evts)-1
   date = []
   board = np.zeros(shape=(np.size(evts)),dtype = np.int32)
   VPower = np.zeros(shape=(nevts,6))
   Th = np.zeros(shape=(nevts,3,2))
   Temp = np.zeros(shape=(np.size(evts)))
   maxCoarse = np.zeros(shape=(np.size(evts)))

   TrigRate = np.zeros(shape=(nevts,7))
   j = 0
   for i in range(j,nevts):
      evt = evts[i]
      if len(evt)>60:  #59 is ACK size
           #print 'SLC word'

           evtsplit = evt.split('\n')
           date.append(evtsplit[1])
           IP = evtsplit[2].split(':')[1]
           board[j] = int(IP[-2:]);

	   for k in range(0,6):
             VPower[j,k]=evtsplit[k+3].split(':')[1]

	   for k in range(0,3):
             thsplit = evtsplit[k+9].split(':')
             Th[j,k,0] = thsplit[1].split(' ')[0]
             Th[j,k,1] = thsplit[1].split(' ')[1]

           Temp[j] = evtsplit[12].split(':')[1]

           for k in range(0,7):
             TrigRate[j,k] = evtsplit[13+k].split(':')[1]

	   maxCoarse[j] = evtsplit[20].split(':')[1]

           j = j+1
   nevts = j
   date = np.array(date)
   #d = mx.DateTime.DateTimeFrom(date[0])  # Better parser than standard datetime module
   #d2 = mx.DateTime.DateTimeFrom(date[10])
   #print d.gmtime(),d2.gmtime(),d2-d

   hraw = [hex(int(a)) for a in Temp]  # Transfer back to hexadecimal
   braw = [bin(int(a)) for a in Temp]
   draw = [twos_comp(int(a,16), 13) for a in hraw] #2s complements

   lab = ['Total','Ch1+','Ch2+','Ch3+','Ch1-','Ch2-','Ch3-']
   boards = set(board[np.where(board>0)])
   print 'Boards in run:',list(boards)
   print 'Run start:', date[0]
   for id in boards:
   	   sel = np.where(board == id)
	   date_end = date[sel[0][-1]]
   	   print 'Run stop:',date_end,'for board',id,' (',np.size(sel),'measurements)'

	   pl.figure(1)
   	   pl.plot(Temp[sel])
	   pl.grid(True)
	   pl.xlabel('SLC measurment nb')
	   pl.ylabel('Temperature')
	   pl.title('Temperature')

	   pl.figure(2)
	   for i in range(np.shape(VPower)[1]):
	        sub=321+i
	   	pl.subplot(sub)
		pl.plot(VPower[sel,i][0])
		pl.grid(True)
		if i>3:
  		  pl.xlabel('SLC measurment nb')
		pl.title('Voltage'+str(i+1))

	   pl.figure(3)  #Trig Rate
	   for i in range(np.shape(TrigRate)[1]):
	        if len(boards)>1:  # Subplot if 2 boards
			sub=211+int(id==list(boards)[-1])
			pl.subplot(sub)
		if np.sum(TrigRate[sel,i][0])>0:
	          pl.plot(TrigRate[sel,i][0],label=lab[i])
		#pl.yscale('log')
		pl.grid(True)
		pl.ylabel('TrigRate [Hz]')
		if (len(boards)<2) | (sub==212 ):
		        pl.xlabel('SLC measurment nb')
   		pl.legend()
		pl.title('Board '+str(id))

	   pl.figure(4)  #MaxCoarse
	   pl.subplot(211)
	   pl.plot(maxCoarse[sel])
	   pl.subplot(212)
	   pl.hist(maxCoarse[sel])
	   print 'maxCoarse counter = ',np.mean(maxCoarse[sel]),'+-',np.std(maxCoarse[sel]),' vs 1.25e+8'

   pl.draw()
   pl.show()

   return

def twos_comp(val, bits):
   """compute the 2's compliment of int value val"""
   if (val & (1 << (bits - 1))) != 0: # if sign bit is set e.g., 8bit: 128-255
     val = val - (1 << bits)	    # compute negative value
   return val



if __name__ == '__main__':
     loop(sys.argv[1])
