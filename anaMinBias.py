# Script to reduce minBias data
# and plot resulting file minBias_b[ID].txt

import os
import time
import sys
import math
from os.path import expanduser

import ephem
import datetime
import numpy as np
import matplotlib
import pylab as pl
pl.ion()

def loopRuns(boardID,startrun,endrun):
  print "Calling loopRuns(). Will analyse minBias for board {0} between R{1} and R{2}.".format(boardID,startrun,endrun)
  time.sleep(1)
  
  datadir = "/home/pastsoft/data"
  #  datadir = "/home/martineau/GRAND/GRANDproto35/data/ulastai/"
  for run in range(int(startrun),int(endrun)+1):
    # Build run name
    filename = datadir+'M'+str(run)+'_b05.data.txt'
    loopEvents(filename,boardID)

def loopEvents(filename,boardID):

   timeslice = 60  # Time duration for the data integration (mins)
   
   resfile = 'minBias_b'+boardID+'.txt'  # Output file
   reso = open(resfile,'ab')
   a = np.loadtxt(resfile)
   try:
     tfmax =  a[-1,0]   
   except IndexError:  # When file is empty 
     tfmax = 0


   if os.path.isfile(filename) is False:
     print 'File ',filename,'does not exist. Abort'
     return
   print 'Scanning minBias datafile',filename
   with open(filename,"r") as f:
   	   evts = f.read().split('-----------------')

   nevts = len(evts)-1
   print 'Number of events:',nevts  # First element does not count
   time.sleep(1)
   
   # Arrays initialization
   unixsecs = []
   data = list()
   
   # Loop on events
   j = 0;  # Index of array filling (because date & data are "append")
   for i in range(0,nevts-1):
           if float(i)/100 == int(i/100):
	   	print 'Fetching event',i,'/',nevts
           evt = evts[i+1]
   	   evtsplit = evt.split('\n')
 
   	   if np.size(evtsplit)>8:   # Event is of normal size
		   
		   IP = evtsplit[2][3:]
		   board = int(IP[-2:]);
		   if board != int(boardID):
		     #print 'This is board {0}, skiping it (analysing board {1} only)'.format(board,boardID)
		     continue
		   
		   # Now reducing data
		   # Time info
		   date = evtsplit[1]
		   thisDatetime = datetime.datetime.strptime(date, '%a %b %d %H:%M:%S %Y GMT')  # Build DateTime object
      		   utcsec = time.mktime(thisDatetime.timetuple())  # Build UTC second
		   if utcsec<=tfmax: # Only looking at data more recent than already present in minBias_b[ID].txt
		     print 'Older data than in {0}, skiping it.'.format(resfile)
		     continue		   
                   unixsecs.append(utcsec)
		                      
		   # Data
   		   raw=evtsplit[9:][:]  #raw data
   		   raw2 = raw[0].split(" ") # Cut raw data list into samples
		   raw2 = raw2[0:np.size(raw2)-1]   # Remove last element (empty)
		   hraw2 = [hex(int(a)) for a in raw2]  # Transfer back to hexadecimal
   		   draw = [twos_comp(int(a,16), 12) for a in hraw2] #2s complements
		   nsamples = len(draw)/4  # Separate data to each channel
		   thisEvent = np.reshape(draw,(4,nsamples));
   		   data.append(thisEvent) # Write to data list
		   
		   j = j+1
   	   else:
   		   print 'Error! Empty event',i

   #
   nevtsb = np.shape(data)[0]
   sig = np.asarray([])
   tt = np.asarray([])
   lst = np.asarray([])
   # Local sideral time info @ Ulastai
   ulastai = ephem.Observer();
   ulastai.long = ephem.degrees("86.71")
   ulastai.lat = ephem.degrees("42.95")
   ulastai.elevation = 2650;
       
   # Analysing 
   if nevtsb==0:  # No data read for this antenna in this run     
     return
     
   print '{0} new data points in {1}.\nNow adding them to minBias result file {2}.'.format(nevtsb,filename,resfile)
   
   # Now group events
   comEventX = []
   comEventY = []
   comEventZ = []
   j = 0
   i = 1    
   istart = 0
   
   while i < min(nevtsb,500000):  #Loop on all events from this board
     if float(i)/100 == int(i/100):
       print 'Analysing event',i,'/',nevtsb
     ind = i
     #print i,nevtsb,np.shape(data)
     thisEvent = data[ind]  
     if unixsecs[i]-unixsecs[istart]<timeslice:  # Close event in time: merge them
       #print 'Adding event',i
       comEventX = np.concatenate((comEventX,thisEvent[0][:]),axis=0)
       comEventY = np.concatenate((comEventY,thisEvent[1][:]),axis=0)
       comEventZ = np.concatenate((comEventZ,thisEvent[2][:]),axis=0)
     else:  # Time slice now ended. 
       # Saving last slice to sig array
       ulastai.date = datetime.datetime.fromtimestamp(unixsecs[istart])
       sig = np.append(sig,[unixsecs[i],ulastai.sidereal_time(),np.mean(comEventX), np.mean(comEventY), np.mean(comEventZ)])   
       
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
       
       # Reset parameters for next slice
       istart = i  
       j = j+1
       comEventX = []
       comEventY = []
       comEventZ = []
       
     i = i+1;
   
   # Add last event     
   if j == 0:
     ulastai.date = datetime.datetime.fromtimestamp(unixsecs[istart])    
     sig = np.append(sig,[unixsecs[-1],ulastai.sidereal_time(),np.mean(comEventX), np.mean(comEventY), np.mean(comEventZ)])	 

   # End of loop on all data. Now write to file.
   sig = sig.reshape(np.size(sig)/5,5)
   np.savetxt(reso, sig)


def displayGalVar(boardID):
   home = expanduser("~")
   resdir = home+"/GRAND/GRANDproto35/data/ulastai/"
   resfile = resdir+"minBias_b"+str(boardID)+".txt"
   print "Calling displayGalVar(). Will display minBias result file {0}".format(resfile)
   
   sd,sm,sy=28,9,2018  # Start day,month,year
   ed,em,ey=1,1,2024  # End day,month,year
   print "Period displayed: {0}/{1}/{2}-{3}/{4}/{5}".format(sd,sm,sy,ed,em,ey)
   startwindow=(datetime.datetime(sy,sm,sd)-datetime.datetime(1970,1,1)).total_seconds()
   endwindow=(datetime.datetime(ey,em,ed)-datetime.datetime(1970,1,1)).total_seconds()
   time.sleep(1)
   
   # Load data
   adcF = 2./4096 # 12bits 2V 
   reso = open(resfile,'ab')
   a = np.loadtxt(resfile)
   t = a[:,0]
   lst = a[:,1]
   x = a[:,2]*adcF
   y = a[:,3]*adcF
   z = a[:,4]*adcF
   sel = np.where((t<endwindow) & (t> startwindow))
   t = t[sel]
   lst = lst[sel]
   x = x[sel]
   y = y[sel]
   z = z[sel] 
   res = np.asarray([x,y,z]).reshape(3,np.size(x))
   
   nticks = 8
   ind = np.linspace(min(t),max(t),nticks)
   date = [datetime.datetime.fromtimestamp(ux).strftime('%m/%d') for ux in ind]
   lsth = lst*24./(2*np.pi); # Now in 0:24h range

  
   for k in [0, 1, 2]:
     thisSig = res[k,:]
     
     pl.figure(1)
     pl.plot(t,thisSig-0.05*k,'+')
     pl.xticks(np.linspace(min(t),max(t),nticks), date)
     pl.xlim(min(t)-1,max(t)+1)
     pl.xlabel('Date [Month/Day]',size='large')
     pl.ylabel('Mean ADC voltage [V]',size='large')
     pl.grid(True)
     pl.show()
   
     pl.figure(2)  # Folded
     subp = 311+k
     pl.subplot(subp)
     pl.plot(lsth,thisSig,'+')
     pl.plot(lsth[-1],thisSig[-1],'o')
     pl.xticks([h for h in range(0,24,2)], ['%i:00'%h for h in range(0,24,2)])
     pl.xlim(0,24)
     pl.title('Channel {0} - {1}'.format(k,resfile))
     pl.ylabel('Mean baseline (V)')
     if k == 2:
       pl.xlabel('Local Sideral Time (h)')
     pl.grid(True)
     pl.show()
   pl.savefig('{0}.png'.format(resfile))
  
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
     
     #loopRuns(sys.argv[1],sys.argv[2],sys.argv[3])
     displayGalVar(sys.argv[1])
     
