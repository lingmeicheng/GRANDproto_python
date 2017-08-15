import os
import time
import sys
import math

import ephem
import datetime
import numpy as np
import matplotlib
#matplotlib.use('Agg')
import pylab as pl


def loopEvents(filename):

   boardID = 18
   timeslice = 600
   
   resfile = 'minBias.txt'
   reso = open(resfile,'ab')
   a = np.loadtxt(resfile)
   if np.size(a)>0:
     tfmax =  a[-1,0]
   else:
     tfmax = 0
   
   datafile = '../data/'+filename
   print 'Scanning minBias datafile',datafile
   with open(datafile,"r") as f:
   	   evts = f.read().split('-----------------')

   nevts = len(evts)-1
   print 'Number of events:',nevts  # First element does not count
   time.sleep(1)
   unixsecs = []
   data = list()
  
   #pl.ion()

   j = 0;  # Index of array filling (because date & data are "append")
   for i in range(0,nevts-1):
           if float(i)/100 == int(i/100):
	   	print 'Fetching event',i,'/',nevts
           evt = evts[i+1]
   	   evtsplit = evt.split('\n')
 
   	   if np.size(evtsplit)>8:   # Event is of normal size
		   
		   IP = evtsplit[2][3:]
		   board = int(IP[-2:]);
		   if board != boardID:
		     print 'Board {0}, skip it (analysing board {1} only'.format(board,boardID)
		     continue
		   
		   date = evtsplit[1]
		   a = evtsplit[1].split()
      		   year = int(a[4])
      		   if a[1]=='May':
		     month = 5
		   elif  a[1]=='June':   
		     month = 6
      		   day = int(a[2])
      		   h = int(a[3].split(':')[0])
      		   mn = int(a[3].split(':')[1])
      		   s = int(a[3].split(':')[2])
      		   dstrf = '{0}/{1}/{2} {3}:{4}:{5}'.format(year,month,day,h,mn,s)
      		   #print 'dstrf:',dstrf
      		   d = datetime.datetime(year,month,day,h,mn,s)
      		   sec = time.mktime(d.timetuple())
		   if sec<=tfmax:
		     #print 'Date already in {0}, skip data.'.format(resfile)
		     continue
		
		   # Now only fresh data on valid board
		   unixsecs.append(sec)                   
		   # Data
   		   raw=evtsplit[9:][:]  #raw data
   		   raw2 = raw[0].split(" ") # Cut raw data list into samples
		   raw2 = raw2[0:np.size(raw2)-1]   # Remove last element (empty)
		   hraw2 = [hex(int(a)) for a in raw2]  # TRansfer back to hexadecimal
   		   draw = [twos_comp(int(a,16), 12) for a in hraw2] #2s complements
		   nsamples = len(draw)/4  # Separate data to each channel
		   thisEvent = np.reshape(draw,(4,nsamples));
   		   data.append(thisEvent) # Write to data list
		   
		   j = j+1
   	   else:
   		   print 'Error! Empty event',i

   nevtsb = np.shape(data)[0]
   sig = np.asarray([])
   tt = np.asarray([])
   lst = np.asarray([])
   # Local sideral time info
   ulastai = ephem.Observer();
   ulastai.long = ephem.degrees("86.71")
   ulastai.lat = ephem.degrees("42.95")
   ulastai.elevation = 2650;

   for b in [0]:
     
     # Group events
     comEventX = []
     comEventY = []
     comEventZ = []
     
     print '{0} new data points in {1}.\nNow adding them to minBias result file {2}.'.format(nevtsb,datafile,resfile)
     
     j = 0
     i = 1    
     istart = 0
     while i < min(nevtsb,500000):  #Loop on all events from this board
       if float(i)/100 == int(i/100):
	 print 'Analysing event',i,'/',nevtsb
       ind = i
       #print i,nevtsb,np.shape(data)
       thisEvent = data[ind]  
       if unixsecs[i]-unixsecs[istart]<timeslice:  # Close event in time: merge
         #print 'Adding event',i
         comEventX = np.concatenate((comEventX,thisEvent[0][:]),axis=0)
         comEventY = np.concatenate((comEventY,thisEvent[1][:]),axis=0)
         comEventZ = np.concatenate((comEventZ,thisEvent[2][:]),axis=0)
       else:
         #print 'Done!'
	 #print b,istart,i,unixsecs[istart],unixsecs[i],np.mean(comEventX)
	 ulastai.date = datetime.datetime.fromtimestamp(unixsecs[istart])
	 sig = np.append(sig,[unixsecs[i],ulastai.sidereal_time(),np.mean(comEventX), np.mean(comEventY), np.mean(comEventZ)])	 
	 #lst = np.append(lst,ulastai.sidereal_time())  # LST is an PyEphem angle (in hours)
	 #print ulastai.date,ephem.hours(lst[j])

	 
	 if 0:  #DISPLAY
 	   pl.figure(1)
 	   pl.subplot(311)
 	   pl.plot(comEventX)
 	   pl.subplot(312)
 	   pl.plot(comEventY)
 	   pl.subplot(313)
 	   pl.plot(comEventZ)
 	   #pl.show()
 	   raw_input()
 	   pl.close(1)
	 
	 istart = i  
	 j = j+1
	 comEventX = []
         comEventY = []
         comEventZ = []
	 
       i = i+1;
   
   # End of loop on all data. Now write to file.
   sig = sig.reshape(np.size(sig)/5,5)
   np.savetxt(reso, sig)

def displayGalVar():
   pl.ion()
   font = {'family' : 'normal','weight' : 'bold','size'   : 22}
 
   pl.rc('font', **font)
   
   adcF = 2./4096 # 12bits 2V 
   resfile = 'minBias.txt'
   reso = open(resfile,'ab')
   a = np.loadtxt(resfile)
   t = a[:,0]
   lst = a[:,1]
   x = a[:,2]*adcF
   y = a[:,3]*adcF
   z = a[:,4]*adcF
   sel = np.where(t<1496000000)
   t = t[sel]
   lst = lst[sel]
   x = x[sel]
   y = y[sel]
   z = z[sel] 
   res = np.asarray([x,y,z]).reshape(3,np.size(x))
   
   nticks = 8
   ind = np.linspace(min(t),max(t),nticks)
   date = [datetime.datetime.fromtimestamp(ux).strftime('%Y/%m/%d') for ux in ind]
   lsth = lst*24./(2*np.pi); # Now in 0:24h range

  
   for k in [0, 1, 2]:
     thisSig = res[k,:]
     
     pl.figure(1)
     pl.plot(t,thisSig-0.05*k,'+')
     pl.xticks(np.linspace(min(t),max(t),nticks), date)
     pl.xlabel('Date [Year/Month/Day]',size='large')
     pl.ylabel('Mean ADC voltage [V]',size='large')
     pl.grid(True)
     pl.show()
   
     pl.figure(2)  # Folded
     subp = 311+k
     pl.subplot(subp)
     pl.plot(lsth,thisSig,'+')
     pl.xticks([h for h in range(0,24,2)], ['%i:00'%h for h in range(0,24,2)])
     pl.xlim(0,24)
     pl.title('Channel {0} - {1}'.format(k,resfile))
     pl.ylabel('Mean bline [V]')
     if k == 2:
       pl.xlabel('Local Sideral Time [h]')
     pl.grid(True)
     pl.show()
     pl.savefig('minbias.png')
  
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
     #loopEvents(sys.argv[1])
     displayGalVar()
     pl.ion()
     