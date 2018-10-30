# Analyses calib events from tests
# Taken from anaData.py
import os
import time
import sys
import math

import numpy as np
import pylab as pl
from scipy.optimize import curve_fit

def loopEvents(RUNID,folder,BOARDID,TYPE):
   pl.ion()
   DISPLAY = 1
   nch = 4   #Nb of channels
   print 'DISPLAY = ',DISPLAY 
   
   pl.ion()
   if TYPE=='1':
     pre='C'
   elif TYPE=='4':
     pre='M'
   elif TYPE=='2':
     pre='R'
   else:
     print 'Wrong TYPE!'
     
   datafile = folder+'/'+pre+str(RUNID)+'_b'+str(BOARDID)+'.data'  
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
   imax = np.zeros(shape=(nevts,nch),dtype=int)
   Amax = np.zeros(shape=(nevts,nch))
   mub = np.zeros(shape=(nevts,nch))
   sigb = np.zeros(shape=(nevts,nch))
   
   data = list()

   j = 0;  # Index of array filling (because date & data are "append")
   for i in range(1,nevts+1):  
   	   if float(i)/100 == int(i/100):
	   	print 'Event',i,'/',nevts
   	   evt = evts[i]
   	   evtsplit = evt.split('\n')
	   #print evt
   	   if np.size(evtsplit)>8:   # Event is of normal size
	           #print evtsplit[3]
   		   date.append(evtsplit[1])
		   IP = evtsplit[2][3:]
		   board[j] = int(IP[-2:]);
		   		
   		   TS2[j]=int(evtsplit[3][4:])  # time elapsed since last PPS (125MHz clock <=> 8ns counter)
   		   #print TS2
		   tt=int(evtsplit[4][11:])  # phase in 8ns slot fr trigger
   		   TS1Trig[i] = tt
		   #TS1Trig[i] = get_1stone(hex(tt))
		   tpps=int(evtsplit[5][7:]) 
		   #TS1PPS[j]=get_1stone(hex(tpps))  # phase in 8ns slot for PPS
   		   TS1PPS[j]=tpps
		   SSS[j]=int(evtsplit[6][4:])  # Elapsed seconds since start
   		   EvtId[j] = int(evtsplit[7][3:])
   		   TrigPattern[j] = int(evtsplit[8][12:])
   		   # Data
   		   raw=evtsplit[9:][:]  #raw data
   		   raw2 = raw[0].split(" ") # Cut raw data list into samples
		   raw2 = raw2[0:np.size(raw2)-1]   # Remove last element (empty)
	           hraw2 = [hex(int(a)) for a in raw2]  # Transfer back to hexadecimal
		   draw = [twos_comp(int(a,16), 12) for a in hraw2] #2s complements		   
                   draw = np.array(draw)*1./2048  # in Volts
		   
   		   nsamples = len(draw)/4  # Separate data to each channel
   		   offset = int(nsamples/2.0)
   		   #print nsamples,"samples per channel --> offset = ",offset
		   thisEvent = np.reshape(draw,(4,nsamples));
   		   data.append(thisEvent) # Write to data list
		   if DISPLAY:
		     print 'Event ',j, 'at date',date[j]
		     t = np.array(range(np.shape(thisEvent)[1]))
 		     t = t* 10e-3  #in mus
 		     pl.figure(j)
		     pl.subplot(221)
 		     pl.plot(t[3:],thisEvent[0][3:])
 		     pl.xlabel('Time ($\mu$s)')
 		     pl.ylabel('X amplitude (V)')
 		     pl.grid(True)
		     pl.subplot(222)
 		     pl.xlabel('Time ($\mu$s)')
 		     pl.ylabel('Y amplitude (V)')
 		     pl.plot(t[3:],thisEvent[1][3:])
 		     pl.grid(True)
		     pl.subplot(223)
 		     pl.plot(t[3:],thisEvent[2][3:])
 		     pl.xlabel('Time ($\mu$s)')
 		     pl.ylabel('Z amplitude (V)')
 		     pl.grid(True)
		     pl.subplot(224)
 		     pl.plot(t[3:],thisEvent[3][3:])
 		     pl.xlabel('Time ($\mu$s)')
 		     pl.ylabel('Calibrator amplitude (V)')
 		     pl.grid(True)

 		     pl.suptitle('Board {0} Event {1}'.format(board[j],EvtId[j]))
		     
 		     pl.plot(t[3:],thisEvent[3][3:],'s')
  		     xr = t[3:]  #mus
 		     w = 2*np.pi*66.666666  #rad/mus
		     yr = thisEvent[3][3:]
		     fitfunc = lambda xr, a, b, c: a*np.sin(w*xr+b)+c	# Create fit function
		     abeg = float(np.max(yr)-np.min(yr))
		     p, pcov = curve_fit(fitfunc,xr,yr,p0 = [abeg,0.0,0.0])  #Perform fit
 		     print 'Fit results:',p,np.sqrt(np.diag(pcov))
		     xf=np.linspace(xr[0],xr[-1],10000)  # Display fit result wuith nice thinning
		     pl.plot(xf,fitfunc(xf,p[0],p[1],p[2]))
		     
		     pl.show()
		     raw_input()
 		     pl.close(j)
		   
		   for k in range(nch):
		     imax[j,k] = np.argmax(thisEvent[k][3:])+3;  # Skip 1st 3 points because could be left overs from previous events
		     #print k,np.argmax(thisEvent[k][3:]),thisEvent[k][imax[j,k]];
		     #print thisEvent[k][:]
		     #raw_input()
		     Amax[j,k] = thisEvent[k][imax[j,k]];
		     mub[j,k] = np.mean(thisEvent[k][1:offset-5])
		     sigb[j,k] = np.std(thisEvent[k][1:offset-5])
                  

		   j = j+1
   	   else:
   		   print 'Error! Empty event',i
   if TYPE=='1':  
     filename = 'summary_b{}.txt'.format(BOARDID)
     filename = folder + '/' + filename
     fil = open(filename,'a')
     fil.write('*** Calibration: analysis of run {0} \n'.format(RUNID))

     for k in range(nch):
      if k<3:
     	outtxt = 'Channel {0} : bline mean = {1} V, std dev = {2} V \n'.format(k,np.mean((mub[:,k][0])),np.mean((sigb[:,k][0])))
     	print outtxt
     	fil.write(outtxt)
     	#np.savetxt(filename,outtxt,fmt='%s')
      else:
     	outtxt =  'Channel {0} (calibrator): Peak mean = {1} V \n'.format(k,np.mean((Amax[:,k][0])))
     	print outtxt
     	fil.write(outtxt)
     fil.close()
   elif TYPE=='4':
     filename = 'freq_b{}.txt'.format(BOARDID)
     filename = folder + '/' + filename
     fil = open(filename,'a')
     fil.write('### Frequency calibration: analysis of run {0} \n'.format(RUNID))
     for k in range(3):
     	outtxt = 'Channel {0} : bline mean = {1} V, std dev = {2} V \n'.format(k,np.mean((mub[:,k][0])),np.mean((sigb[:,k][0])))
     	print outtxt
     	fil.write(outtxt)
     	#np.savetxt(filename,outtxt,fmt='%s')
     
   fil.close()

     
   j = j+1
  
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
     #print sys.argv
     loopEvents(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
