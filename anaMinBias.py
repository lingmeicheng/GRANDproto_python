import os
import time
import sys
import math

import ephem
import datetime
import numpy as np
import pylab as pl

def loopEvents(filename):

   ulastai = ephem.Observer();
   ulastai.long = ephem.degrees("86.71")
   ulastai.lat = ephem.degrees("42.95")
   ulastai.elevation = 2650;

   datafile = '../data/'+filename
   print 'Scanning',datafile

   with open(datafile,"r") as f:
   	   evts = f.read().split('-----------------')

   nevts = len(evts)-1
   print 'Number of events:',nevts
   time.sleep(1)
   date = []
   board = np.zeros(shape=(np.size(evts)),dtype = np.int32)
   TS2 = np.zeros(shape=(np.size(evts)))
   TS1PPS = np.zeros(shape=(np.size(evts)))
   TS1Trig = np.zeros(shape=(np.size(evts)))
   SSS = np.zeros(shape=(np.size(evts)),dtype = np.int32)
   EvtId = np.zeros(shape=(np.size(evts)),dtype = np.int32)
   TrigPattern = np.zeros(shape=(np.size(evts)))
   imax = np.zeros(shape=(nevts,3))
   Amax = np.zeros(shape=(nevts,3))
   mub = np.zeros(shape=(nevts,3))
   sigb = np.zeros(shape=(nevts,3))
   data = list()
  
   pl.ion()

   j = 0;  # Index of array filling (because date & data are "append")
   for i in range(j,nevts):  
   	   if float(i)/100 == int(i/100):
	   	print 'Event',i,'/',nevts
   	   evt = evts[i]
   	   evtsplit = evt.split('\n')
 
   	   if np.size(evtsplit)>8:   # Event is of normal size
   		   date.append(evtsplit[1])
		   IP = evtsplit[2][3:]
		   board[j] = int(IP[-2:]);
		   		
   		   TS2[j]=int(evtsplit[3][4:])  # time elapsed since last PPS (125MHz clock <=> 8ns counter)
   		   tt=int(evtsplit[4][11:])  # phase in 8ns slot fr trigger
   		   TS1Trig[i] = get_1stone(hex(tt))
		   tpps=int(evtsplit[5][7:]) 
		   TS1PPS[j]=get_1stone(hex(tpps))  # phase in 8ns slot for PPS
   		   SSS[j]=int(evtsplit[6][4:])  # Elapsed seconds since start
   		   EvtId[j] = int(evtsplit[7][3:])
   		   TrigPattern[j] = int(evtsplit[8][12:])
   		   # Data
   		   raw=evtsplit[9:][:]  #raw data
   		   raw2 = raw[0].split(" ") # Cut raw data list into samples
		   raw2 = raw2[0:np.size(raw2)-1]   # Remove last element (empty)
		   hraw2 = [hex(int(a)) for a in raw2]  # TRansfer back to hexadecimal
   		   draw = [twos_comp(int(a,16), 12) for a in hraw2] #2s complements
		   
 
   		   nsamples = len(draw)/4  # Separate data to each channel
   		   offset = nsamples/2.0
   		   #print nsamples,"samples per channel --> offset = ",offset
		   thisEvent = np.reshape(draw,(4,nsamples));
   		   data.append(thisEvent) # Write to data list
		   
		   j = j+1
   	   else:
   		   print 'Error! Empty event',i

   trigtime = np.zeros(shape=(np.size(evts)))
   tdeb = min(SSS[np.where(SSS>0)])
   tend = max(SSS)
   dur = tend-tdeb+1 # Run duration [seconds]
   t = range(dur)
   DataRate = np.zeros(shape=(dur,2))
   TrigRate = np.zeros(shape=(dur,2))
   boards = set(board[np.where(board>0)])
   print 'Boards in run:',boards
   boards = list(boards)
   print 'Run start:', date[0]
   
   sig = np.zeros(shape=(2,3,nevts))
   tt = np.zeros(shape=(1,nevts))
   lst = np.zeros(shape=(nevts),dtype=float)
   
   for b in range(1):
   #range(len(boards)):
     sel = np.where(board == boards[b])
     date_end = date[sel[0][-1]]
     print 'Run stop:',date_end,'for board',boards[b],' (',np.size(sel),'measurements)'
     cor=125.0e6/float(max(TS2[sel]))  # Ratio of expected nb of counts in 1s to actually measured  => correction to clock
     print 'Correction factor for 125MHz clock for board',boards[b],':',cor
     # Build trig time
     trigtime[sel] = SSS[sel]+(TS2[sel]*4+TS1PPS[sel]-TS1Trig[sel])*2e-9*cor  #second. Use same SSS for both cards
     
     # Group events
     comEventX = []
     comEventY = []
     comEventZ = []

     sel = sel[:][0]
     bSSS = SSS[sel]  #SSS for this board
     
     print 'Computing minBias curve for board',boards[b]
     j = 0
     i = 1    
     istart = 0
     while i <len(sel):  #Loop on all events from this board
       #print istart,i,j
       ind = sel[i]
       thisEvent = data[ind]  
       #print bSSS[i],bSSS[i-1],bSSS[i]-bSSS[i-1]
       if bSSS[i]-bSSS[istart]<1200:  # Close event in time: merge
         #print 'Adding event',i
         comEventX = np.concatenate((comEventX,thisEvent[0][:]),axis=0)
         comEventY = np.concatenate((comEventY,thisEvent[1][:]),axis=0)
         comEventZ = np.concatenate((comEventZ,thisEvent[2][:]),axis=0)
       else:
         #print 'Done!'
	 #print b,i,np.mean(comEventX)
	 sig[b,0,j] = np.mean(comEventX)
	 sig[b,1,j] = np.mean(comEventY)
	 sig[b,2,j] = np.mean(comEventZ)
	 #tt[j] = bSSS[istart]
	 
	 dstr = date[sel[istart]]
	 year = 2016
         month=12
	 day=dstr[8:10]
	 h = dstr[11:13]
	 mn = dstr[14:16]
	 s = dstr[17:19]	 
	 dstrf = '2016/12/{0} {1}:{2}:{3}'.format(day,h,mn,s)
	 ulastai = ephem.Observer();
	 ulastai.long = ephem.degrees("86.71")
         ulastai.lat = ephem.degrees("42.95")
         ulastai.elevation = 2650;
	 ulastai.date = ephem.Date(dstrf)
	 lst[j] = ulastai.sidereal_time()  # LST is an PyEphem angle (in hours)
	 print dstrf,ulastai.sidereal_time(),ephem.hours(lst[j])
	 
	 
	 if 0:  #DISPLAY
 	   pl.figure(1)
 	   pl.subplot(311)
 	   pl.plot(comEventX)
 	   pl.subplot(312)
 	   pl.plot(comEventY)
 	   pl.subplot(313)
 	   pl.plot(comEventZ)
 	   pl.show()
 	   raw_input()
 	   pl.close(1)
	 
	 istart = i  
	 j = j+1
	 comEventX = []
         comEventY = []
         comEventZ = []
	 
       i = i+1;
   
   # Display
   for b in range(len(boards)):
     pl.figure(2+b)
     for k in [0, 1, 2]:
       thisSig = sig[b,k,:]
       sigclean = thisSig[np.where(thisSig<0)]
       lstclean=lst[np.where(thisSig<0)]
       tt = []
       for ii in range(0,len(lstclean),1):
	 tt.append(str(ephem.hours(lstclean[ii]))[0:5])       
       subp = 311+k
       pl.subplot(subp)
       pl.plot(sigclean)
       pl.xticks(range(len(tt)),tt,size='small')
       pl.title('Board {0} - Channel {1} - {2}'.format(boards[b],k,filename))
       pl.ylabel('Mean bline [LSB]')
       if k == 2:
         pl.xlabel('Local Sideral Time [h]')
       pl.grid(True)
   pl.show()
 
  
def get_1stone(val):
    if val == '0x1':
    	return 0
    if val == '0x3':
    	return 1
    if val == '0x7':
    	return 2
    if val == '0xf':
    	return 3
    if val == '0x1f':
    	return 4
    if val == '0x3f':
    	return 5
    if val == '0x7f':
    	return 6
    if val == '0xff':
    	return 7
    return 8	
	
	
    
def twos_comp(val, bits):
    """compute the 2's compliment of int value val"""
    if (val & (1 << (bits - 1))) != 0: # if sign bit is set e.g., 8bit: 128-255
        val = val - (1 << bits)        # compute negative value
    return val   
    
   

if __name__ == '__main__':
     loopEvents(sys.argv[1])
