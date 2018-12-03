# Base script to analyse standard dat
# OMH AUg 29, 2018

import os
import time
import sys
import math
import mx.DateTime
import numpy as np
import pylab as pl
from scipy.optimize import curve_fit

# Binary dat ainterpreter: pyef
# module. Do ./setup.py build & ./setup.py install if not loadable [email Gu Junhua Sept 27 2017]
#import pyef


DISPLAY = 0
print 'DISPLAY = ',DISPLAY 
pl.ion()

def loopEvents(filename,SLW_filename):
   if os.path.isfile(filename) is False:
     print 'File ',filename,'does not exist. Aborting.'
     return
  
   if os.path.isfile(SLW_filename) is False:
     print 'File ',filename,'does not exist. Aborting.'
     return
  
   datafile = filename.split('/')
   TYPE = datafile[-1][0]
   if TYPE == 'C':
     nch = 4   # Plot calibrator channel
   else:
     nch = 3

   # Read data
   print 'Scanning data file',filename
   with open(filename,"r") as f:
   	   evts = f.read().split('-----------------')
   nevts = len(evts)-1
   print 'Number of events:',nevts
   time.sleep(1)
   
   # Initialize arrays
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
   freq=50e6  # Sampling Frequency
   dt=1.0/freq
   tdeb = 0
   tend = 0

   # Loop on events
   j = 0;  # Index of array filling (because date & data are "append")
   for i in range(1,nevts+1):  
   	   if float(i)/100 == int(i/100):
	   	print 'Event',i,'/',nevts
   	   evt = evts[i]
   	   evtsplit = evt.split('\n')
	   if np.size(evtsplit)>8:   # Event is of normal size
#	           date.append(evtsplit[1])
                   d = mx.DateTime.DateTimeFrom(evtsplit[1])
                   date.append(d)
		   IP = evtsplit[2][3:]
		   board[j] = int(IP[-2:]);
		   		
   		   TS2[j]=int(evtsplit[3][4:])  # time elapsed since last PPS (125MHz clock <=> 8ns counter)
   		   tt=int(evtsplit[4][11:])  # phase in 8ns slot fr trigger
   		   TS1Trig[i] = tt
		   tpps=int(evtsplit[5][7:]) 
		   TS1PPS[j]=tpps
		   SSS[j]=int(evtsplit[6][4:])  # Elapsed seconds since start
   		   EvtId[j] = int(evtsplit[7][3:])
   		   TrigPattern[j] = int(evtsplit[8][12:])
   		   # Data
   		   raw=evtsplit[9:][:]  #raw data
   		   raw2 = raw[0].split(" ") # Cut raw data list into samples
		   raw2 = raw2[0:np.size(raw2)-1]   # Remove last element (empty)
   		   if TYPE == "P":  # Pattern run, keep data value in binary format
			draw = [int(a) for a in raw2] 
		   else:  # Other run types
  		        hraw2 = [hex(int(a)) for a in raw2]  # Transfer back to hexadecimal
			draw = [twos_comp(int(a,16), 12) for a in hraw2] #2s complements		   
                   if TYPE != "P":  # Transfer data to volts 
		     draw = np.array(draw)*1./2048  # in Volts
		   
   		   nsamples = len(draw)/4  # draw corresponds to 4 channels
   		   offset = int(nsamples/2.0)  # Offset position at center of waveform
   		   #print nsamples,"samples per channel --> offset = ",offset
		   thisEvent = np.reshape(draw,(4,nsamples));
   		   data.append(thisEvent) # Write to data list
		   #print 'Sampling frequency=',freq,'MHz, time step=',dt,'s'
		   if DISPLAY:
		     print 'Event ',j, 'at date',date[j]
		     t = dt*np.array(range(np.shape(thisEvent)[1]))
 		     t = t* 1e6  #in mus
 		     pl.figure(j)
		     if TYPE == "C":
   		       pl.subplot(221)
 		     else:
		       pl.subplot(311)
 		     pl.plot(t[3:],thisEvent[0][3:])
 		     pl.xlim(t[3],max(t))
 		     pl.xlabel('Time ($\mu$s)')
 		     if TYPE == "P":
		       pl.ylabel('LSB')
		     else:
		       pl.ylabel('Amplitude (V)')
 		     pl.grid(True)
		     if TYPE == "C":
   		       pl.subplot(222)
 		     else:
		       pl.subplot(312)
 		     pl.xlabel('Time ($\mu$s)')
 		     pl.xlim(t[3],max(t))
 		     if TYPE == "0":
		       pl.ylabel('LSB')
		     else:
		       pl.ylabel('Amplitude (V)')
		     pl.plot(t[3:],thisEvent[1][3:])
 		     pl.grid(True)
		     if TYPE == "C":
   		       pl.subplot(223)
 		     else:
		       pl.subplot(313)
 		     pl.plot(t[3:],thisEvent[2][3:])
                     pl.xlim(t[3],max(t))
 		     pl.xlabel('Time ($\mu$s)')
 		     if TYPE == "P":
		       pl.ylabel('LSB')
		     else:
		       pl.ylabel('Amplitude (V)')
 		     pl.grid(True)
		     if TYPE == "C":
 		       pl.subplot(224)
 		       pl.plot(t[3:],thisEvent[3][3:])
 		       pl.xlabel('Time ($\mu$s)')
 		       pl.ylabel('Amplitude (V)')

 		     pl.grid(True)

 		     pl.suptitle('Board {0} Event {1}'.format(board[j],EvtId[j]))
		     
		     if TYPE == "C":  # Plotting calibrator signal in Calibration mode
 		       pl.plot(t[3:],thisEvent[3][3:],'s')
		       # Fit calibration signal with sine wave
  		       xr = t[3:]  #mus
 		       w = 2*np.pi*66.666666  #rad/mus
		       yr = thisEvent[3][3:]
		       fitfunc = lambda xr, a, b, c: a*np.sin(w*xr+b)+c   # Create fit function
		       abeg = float(np.max(yr)-np.min(yr))
		       p, pcov = curve_fit(fitfunc,xr,yr,p0 = [abeg,0.0,0.0])  #Perform fit
 		       print 'Fit results:',p,np.sqrt(np.diag(pcov))
		       xf=np.linspace(xr[0],xr[-1],10000)  # Display fit result wuith nice thinning
		       pl.plot(xf,fitfunc(xf,p[0],p[1],p[2]))
		     
		     pl.show()
		     raw_input()
 		     pl.close(j)
		   
		   # Assemble stat for summary plots
		   for k in range(nch):
		     imax[j,k] = np.argmax(thisEvent[k][3:])+3;  # Skip 1st 3 points because could be left overs from previous events
		     Amax[j,k] = thisEvent[k][imax[j,k]];
		     mub[j,k] = np.mean(thisEvent[k][1:offset-5])
		     sigb[j,k] = np.std(thisEvent[k][1:offset-5])                  

		   j = j+1
   	   else:
   		   print 'Error! Empty event',i

   if TYPE == "P":
     return
   
   # Now display summary plots
   trigtime = np.zeros(shape=(np.size(evts))) 
   deltat = np.zeros(shape=(np.size(evts)))   #Delta_t to previous event on same antena
   
   timein = np.where(SSS>0)  
   if np.size(timein) > 0:  # GPS timing info is available
     tdeb = min(SSS[timein])
     tend = max(SSS)
   dur = tend-tdeb+1 # Run duration [seconds]
   t = range(dur)
   boards = set(board[np.where(board>0)])
   DataRate = np.zeros(shape=(dur,len(boards)))
   TrigRate = np.zeros(shape=(dur,len(boards)))
   
   print 'Run start:', date[0]
   print 'Boards in run:',list(boards)
   j = 0
   for id in boards:  # Loop on all boards in run
     sel = np.where(board == id)
     date_end = date[sel[0][-1]]
     print 'Run stop:',date_end,'for board',id,' (',np.size(sel),'measurements)'
     if np.size(timein) > 0:
       # To be implemented: read MaxCoarse info from slow control data and define correction factor accordingly cor = 125e6/MaxCoarse
       # See anaTiming.py line 133 and following for guidance

#------------------------------------------
        # Grab MaxCoarse info from SLC run.
        [d,maxCoarseRaw] = getMaxCoarse(SLW_filename)
        
        d = np.array(d)
        maxCoarse = np.zeros((np.size(sel)))
        
        for k in range(np.size(sel)): 
     	  idift = np.argmin(np.abs(date[k]-d))
          maxCoarse[k] = maxCoarseRaw[idift]
        cor=125e6/(maxCoarse+1)
        
#       cor = 1.0
#------------------------------------------------
     else:
       cor=1.0
     print 'Correction factor for 125MHz clock for board',id,':',cor
     # Build trig time
     trigtime[sel] = SSS[sel]+(TS2[sel]*4+TS1PPS[sel]-TS1Trig[sel])*2e-9*cor  #second. 
     
     deltat[sel] = np.diff(trigtime[sel])
     # Compute trig rate
     for i in range(dur):
	ts = tdeb+i
	thisSec = np.where(SSS[sel]==ts)
	thisEvtId = EvtId[sel]
	if np.size(thisSec) > 0:
	  thisSec=thisSec[:][0]
	  DataRate[i,j] = np.size(thisSec)
	  TrigRate[i,j] = thisEvtId[thisSec[-1]]-thisEvtId[thisSec[0]]+1   #Nb of events rigged in that second --> trigrate
     
     pl.figure(2)
     pl.plot(t,TrigRate[:,j],label='Trig rate - Board '+str(id))     
     pl.plot(t,DataRate[:,j],label='Data rate - Board '+str(id))
     pl.grid(True)
     pl.xlabel('Run time (s)')
     pl.ylabel('Data rate (Hz)')
     pl.legend()
     pl.title('Data rate')    
     
     # Displays
     for k in range(nch):
      if 1:
       good = np.where( (imax[sel,k][0]>104) & (imax[sel,k][0]<108))
       abline = np.where( (Amax[sel,k][0]<0))
       azero = np.where( (Amax[sel,k][0]==0))
       print 'Channel',k,': good events=',np.size(good),'/',np.size(sel),'=',float(np.size(good))/np.size(sel)
       print 'Channel',k,': Max at zero=',np.size(azero),'/',np.size(sel),'=',float(np.size(azero))/np.size(sel)
       print 'Channel',k,': Max < zero=',np.size(abline),'/',np.size(sel),'=',float(np.size(abline))/np.size(sel)
       
       pl.figure(id*100+21+k)
       
       pl.subplot(231)
       pl.hist(mub[sel,k][0],offset*2)
       pl.xlabel('Baseline mean')
       pl.title('Board {0}'.format(id))
       pl.grid(True)

       pl.subplot(235)
       pl.plot(mub[sel,k][0],'+')
       pl.plot(Amax[sel,k][0],'o')
       pl.xlabel('Event ID')
       pl.ylabel('Mean amp (bline & max)')
       pl.title('Board {0}'.format(id))
       pl.grid(True)

       pl.subplot(234)
       pl.xlabel('Index of signal max')
       pl.hist(imax[sel,k][0],offset*2)
       pl.title('Board {0}'.format(id))
       pl.grid(True)
       
       pl.subplot(236)
       pl.xlabel('Max amplitude')
       pl.hist(Amax[sel,k][0],offset*2)
       pl.title('Board {0}'.format(id))
       pl.grid(True)
              
       pl.subplot(232)
       diffAmp = Amax[sel,k][0]-mub[sel,k][0]
       pl.hist(sigb[sel,k][0],offset*2)
       pl.xlabel('Bline std dev')
       pl.title('Board {0}'.format(id))
       pl.grid(True)
       
       print 'Channel',k,': bline @ ',np.mean((mub[sel,k][0])),'pm',np.std((mub[sel,k][0])),'V. Std dev=',np.mean((sigb[sel,k][0])),'V'
       print 'Channel',k,': Peak @ ',np.mean((Amax[sel,k][0])),'V, std dev=',np.std((Amax[sel,k][0])),'V, rel error=',np.std((Amax[sel,k][0]))/np.mean((Amax[sel,k][0]))*100,'%'
       print 'Channel',k,': Peak - bline @ ',np.mean((diffAmp)),'V, std dev=',np.std((diffAmp)),'V, rel error=',np.std((diffAmp))/np.mean((diffAmp))*100,'%'
       
       pl.subplot(233)
       pl.plot(mub[sel,k][0],sigb[sel,k][0],'+')
       pl.xlabel('Baseline mean')
       pl.ylabel('Bline std dev')
       pl.title('Board {0}'.format(id))
       pl.grid(True)
       
     j = j+1

   
   sel = np.where(trigtime>0)  #GPS time info present
   if np.size(sel)>0:
     first = trigtime[sel[0][0]] 
     last = trigtime[sel[0][-1]]
     dur = last-first
     print 'Nevents = ',np.size(sel)
     print 'Duration [s]',first,last,dur   
     rate = np.size(sel)/dur
     print 'Rate=',rate,'Hz'
     pl.figure(18)
     pl.plot(trigtime)
     pl.xlabel('Evt ID')
     pl.ylabel('Trig time [s]')
     
     if 0:
       # Now check trig rate vs exected
       consigne = 1; #Expected trig period [s]
       deltat = (deltat-consigne)*1e9;  #  [ns]
       seldeltat = np.where(abs(deltat)<2000)  #1mus difference max
       print 'Nevents with ~',consigne,'s time diff =',np.size(seldeltat)
       print 'Time diff offset to',consigne,'s =',np.mean(deltat[seldeltat]),'ns, std dev=',np.std(deltat[seldeltat]),'ns.'
       pl.figure(19)
       pl.subplot(211)
       pl.plot(trigtime[seldeltat],deltat[seldeltat])
       pl.xlabel('Time [s]')
       pl.ylabel('$\Delta$t [ns]')
       pl.grid(True)
       pl.xlim(min(trigtime[seldeltat]),max(trigtime[seldeltat])+1)
       pl.subplot(212)
       pl.hist(deltat[seldeltat],100)
       pl.xlabel('$\Delta$t [ns]')
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
	

def getMaxCoarse(filename):
# retrieve maxCoarse info from bslow control data
    datafile = filename
    print 'Scanning',datafile
    
    #ef=pyef.read_event_file(datafile)
    #pyef.print_event_file(ef)
    
    try:
      with open(datafile,"r") as f:
    	   evts = f.read().split('-----------------')
    except: 
      print "getMaxCoarse: error!!! Could not fond file",filename,". Return 0 instead."
      return 0,0
      
    nevts=len(evts)-1  # All messages in SLC.txt file
    date = []
    maxCoarse = []
    for i in range(nevts):
       evt = evts[i]
       if len(evt)>60:  #59 is ACK size
            evtsplit = evt.split('\n')
            d = mx.DateTime.DateTimeFrom(evtsplit[1])  # Better parser than standard datetime module

            date.append(d)
            maxCoarse.append(evtsplit[20].split(':')[1])
    return date,maxCoarse

 
    
def twos_comp(val, bits):
    """compute the 2's compliment of int value val"""
    if (val & (1 << (bits - 1))) != 0: # if sign bit is set e.g., 8bit: 128-255
        val = val - (1 << bits)        # compute negative value
    return val   
    
   

if __name__ == '__main__':
     if len(sys.argv)!=3:
       print "Usage: >loopEvents path/filename path/SLW_filename"
     else:  
       loopEvents(sys.argv[1],sys.argv[2])
