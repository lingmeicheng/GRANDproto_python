import os
import time
import sys
import math

import numpy as np
import pylab as pl

def loop(filename):
   pl.ion()

   datafile = '../data/'+filename
   print 'Scanning',datafile

   with open(datafile,"r") as f:
   	   evts = f.read().split('-----------------')

   nevts=len(evts)-1
   date = []
   board = np.zeros(shape=(np.size(evts)),dtype = np.int32)
   VPower = np.zeros(shape=(nevts,6))
   Th = np.zeros(shape=(nevts,3,2))
   Temp = np.zeros(shape=(np.size(evts)))
   TrigRate = np.zeros(shape=(nevts,7))
   j = 0
   for i in range(j,nevts):
      evt = evts[i]
      #print evt
      if len(evt)>59:  #59 is ACK size
           #print 'SLC word'
           
           evtsplit = evt.split('\n')
           IP = evtsplit[2].split(':')[1]
	   if IP == "192.168.1.20":
	        #print "Skip 192.168.1.20"
	   	continue
           board[j] = int(IP[-2:]);
           date.append(evtsplit[1])

	   for k in range(0,6):
	     if k == 1:  # VP2: -3.3V
	        fact = -6.9/2.2
	     elif k == 2:  # VP3: +3.3V
	        fact = 6.9/2.2	
             else:
	       fact = 24.0/2
	     VPower[j,k]=float(evtsplit[k+3].split(':')[1])*5.0/pow(2,12)*fact  
           
	   for k in range(0,3):
             thsplit = evtsplit[k+9].split(':')
             Th[j,k,0] = thsplit[1].split(' ')[0]
             Th[j,k,1] = thsplit[1].split(' ')[1]
	     
           Temp[j] = evtsplit[12].split(':')[1]
	   
           for k in range(0,7):
             TrigRate[j,k] = evtsplit[13+k].split(':')[1]
           j = j+1
	      
   print 'Temp = ', Temp[0:100]
   
   hraw = [hex(int(a)) for a in Temp]  # Transfer back to hexadecimal
   #hraw = []
   #for ii in range(len(Temp)):
   #  hraw.append('1')
   #print hraw  
   draw = [twos_comp(int(a,16), 13)*0.0625 for a in hraw] #2s complements
   print 'Draw=',draw[0:20]
   Temp = np.asarray(draw)
   print 'Temp=',Temp[0:20]
   
   lab = ['Total','Ch1+','Ch2+','Ch3+','Ch1-','Ch2-','Ch3-']   
   boards = set(board[np.where(board>0)])
   print 'Boards in run:',list(boards)
   print 'Run start:', date[0]
   vlab = ['Alim +15V', '-3.3V','+3.3V','LNA1','LNA2','LNA3']
   for id in boards:
   	   sel = np.where(board == id)  
	   date_end = date[sel[0][-1]]
   	   print 'Run stop:',date_end,'for board',id,' (',np.size(sel),'measurements)'
	   
	   pl.figure(1)
   	   pl.plot(Temp[sel[0]])
	   #tt = np.zeros(shape=(1,np.shape(sel)[1]),dtype=str)
	   tt = [];
	   for ii in range(0,np.shape(sel)[1],1):
             #print ii,date[ii][10:16]
	     if float(ii)/20==ii/20:
	       tt.append(str(date[ii][11:16]))
	     else:
	       tt.append("")
	   pl.xticks(range(len(tt)),tt,size='small')      
	   pl.grid(True)
	   pl.xlabel('UT [hh:ss]')
	   pl.ylabel('Board temperature [$^o$C]')
	   pl.title('Board {0}'.format(id))
	   
	   pl.figure(2)
	   for i in range(np.shape(VPower)[1]):
	        sub=321+i
	   	pl.subplot(sub)
		pl.plot(VPower[sel,i][0])
		pl.grid(True)
		if i>3:
  		  pl.xlabel('SLC measurment nb')
		pl.title(vlab[i])
           
		
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
		pl.title('Board '+format(id))
		
   pl.show()
   
   return        
 
def twos_comp(val, bits):
   """compute the 2's compliment of int value val"""
   if (val & (1 << (bits - 1))) != 0: # if sign bit is set e.g., 8bit: 128-255
     val = val - (1 << bits)	    # compute negative value
   return val	
    


if __name__ == '__main__':
     loop(sys.argv[1])
