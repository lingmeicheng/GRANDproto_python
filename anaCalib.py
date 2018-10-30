import os
import time
import sys
import re
import math
import numpy as np
import pylab as pl

def loopEvents(RUNID,boardID,att):
   DISPLAY = 0
   print 'DISPLAY = ',DISPLAY 
   pl.ion()
   filename = "C"+RUNID+"_b"+boardID+".data.txt"
   #datafile = '../data/ulastai/'+filename
   datafile = '/home/martineau/GRAND/GRANDproto35/data/ulastai/'+filename
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
   imax = np.zeros(shape=(nevts,3),dtype=int)
   Amax = np.zeros(shape=(nevts,3))
   mub = np.zeros(shape=(nevts,3))
   sigb = np.zeros(shape=(nevts,3))
   
   j = 0;  # Index of array filling (because date & data are "append")
   for i in range(1,nevts+1):  
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
 		   draw = np.array(draw)*1./2048  # in Volts
   		   nsamples = len(draw)/4  # Separate data to each channel
   		   offset = nsamples/2.0
		   thisEvent = np.reshape(draw,(4,nsamples));
   		   #data.append(thisEvent) # Write to data list ... Not needed here
		   if DISPLAY:
		     print 'Event ',j, 'at date',date[j]
		     t = np.array(range(np.shape(thisEvent)[1]))
 		     t = t* 10e-3  #in mus
 		     pl.figure(j)
 		     pl.subplot(311)
 		     pl.plot(t[3:],thisEvent[0][3:])
 		     pl.ylabel('Amplitude [LSB]')
 		     pl.grid(True)
 		     pl.subplot(312)
 		     pl.ylabel('Amplitude [LSB]')
 		     pl.plot(t[3:],thisEvent[1][3:])
 		     pl.grid(True)
 		     pl.subplot(313)
 		     pl.plot(t[3:],thisEvent[2][3:])
 		     pl.xlabel('Time [mus]')
 		     pl.ylabel('Amplitude [LSB]')
		     
 		     pl.grid(True)
 		     pl.suptitle('Board {0} Event {1}'.format(board[j],EvtId[j]))
 		     pl.show()
 		     raw_input()
 		     pl.close(j)

		   
		   for k in [0,1,2]:
		     nz = np.where(thisEvent[k][:]!=0)
		     imax[j,k] = np.argmax(thisEvent[k][:]);
		     Amax[j,k] = thisEvent[k][imax[j,k]];
		     mub[j,k] = np.mean(thisEvent[k][nz])
		     sigb[j,k] = np.std(thisEvent[k][nz])
		   j = j+1
   	   else:
   		   print 'Error! Empty event',i

   boards = set(board[np.where(board>0)])
   print 'Boards in run:',list(boards)
   j = 0
   m = np.empty([len(boards),3])
   em = np.empty([len(boards),3])
   for id in boards:
     sel = np.where(board == id)     
     for k in [0, 1, 2]:      
       pl.figure(1)
       subpl = 311+k
       pl.subplot(subpl)
       a = mub[sel,k][0]
       pl.hist(a,100)
       if k == 0:
         pl.title('Board {0}'.format(id))

       if k == 2:
       	 pl.xlabel('Mean amplitude')
       pl.grid(True)
       pl.figure(2)
       subpl = 311+k
       pl.subplot(subpl)
       b = sigb[sel,k][0]
       pl.hist(b,100)
       if k == 0:
         pl.title('Board {0}'.format(id))
       if k == 2:
       	 pl.xlabel('Std dev')
       pl.grid(True)
                     
       #Pack up results
       m[j,k] = np.mean(a)
       em[j,k] = np.mean(b)       
       print 'Channel',k,': mean=',m[j,k],'; stddev=',em[j,k]

     j = j+1

    
   return {'m':m, 'em':em}
     

def anaRuns(boardID,runstart,runstop):   # Analyse runs and write result to file
  runs=range(runstart,runstop)  
  
  att = np.zeros(np.shape(runs))
  mm = np.empty([len(runs),3])
  emm = np.empty([len(runs),3])

  for i in range(np.size(runs)):  #loop on runs
    runid = runs[i]
    # Grab attenuation values from config file
    folder="/home/martineau/GRAND/GRANDproto35/data/ulastai"
    file = open(folder+"/C"+str(runid)+"_b"+boardID+".cfg", "r")
    for line in file:
      if re.search("Attr1", line):
    	 att1 = int(line.split()[2])
      if re.search("Attr2", line):
    	 att2 = int(line.split()[2])
    att[i] = att1+att2
    print att1,att2,att[i]
    print att
    
    res = loopEvents(str(runid),boardID,str(att[i]))      
    mm[i,:]=list(res['m'][0])
    emm[i,:]=list(res['em'][0])

  # Write to file
  conc = np.r_[att,mm[:,0],emm[:,0],mm[:,1],emm[:,1],mm[:,2],emm[:,2]]   # Concatenate results
  conc = conc.reshape(7,np.size(att)) # Reshape in line = Att, Measure_Ch[i], Error_Ch[i]
  filename = 'caliboutput_b{0}_R{1}R{2}.txt'.format(boardID,runstart,runstop)
  np.savetxt(filename,conc)  # Write to file

    
def anaRes(boardID,runstart,runstop):    
    
  # Read file
  a = np.loadtxt('caliboutput_b{0}_R{1}R{2}.txt'.format(boardID,runstart,runstop))
  att = a[0,:]
  mm = np.empty([len(att),3])
  emm = np.empty([len(att),3])
  for k in range(3):
    mm[:,k] = a[2*k+1,:]
    emm[:,k] = a[2*k+2,:]
  
  sel = np.where((att>100) & (att<200))  #Fit range
  a = 0.25
  b = -63.5-7  #2*3.3dB attenuator insertion loss
  attdB = att*a+b  # ideal: att = 0dB for att=127 on both attenuators 
  vin = 0.262/2*pow(10,attdB/20) # 10^(attdB/20) attenuation & factor 1/2 because of 3 channels split 
  fig = pl.figure(11)
  for k in range(3):
    pl.subplot(2,1,1)
    pl.plot(att,mm[:,k],'o',markersize=6,label='Channel {0}'.format(k))
    #z = np.polyfit(att[sel],mm[sel,k][0],1)  # Linear fit
    #print 'Channel',k,', slope=',z[0],'LSB/dB Att coef'
    #yth = att*z[0]+z[1]
    #pl.plot(att,yth,'y--')
    pl.subplot(2,1,2)
    pl.plot(att,emm[:,k],'o',markersize=6,label='Channel {0}'.format(k))
  pl.subplot(2,1,1)
  pl.grid(True)
  pl.xlabel('$\Sigma$attenuation indexes')
  pl.ylabel('Mean output level [V]')
  pl.legend(loc='best')
  pl.subplot(2,1,2)
  pl.grid(True)
  pl.xlabel('$\Sigma$attenuation indexes')
  pl.ylabel('StdDev [V]')
  pl.legend(loc='best')

  # Now load external calib data Inout = 66MHz+100mV sine wave + att = 30dB + splitter
  if 0:
    f = np.loadtxt('calibExt.txt')
    runex = f[:,0]
    boardex = f[:,1]
    selex = np.where(boardex==int(boardID))[0]
    Vinex = f[selex,2]*0.001  # mV==>V
    Vindaq = Vinex*pow(10,-30./20)*0.5  # Now apply attenuation (attenuator+splitter)
    
    mex = np.empty([len(selex),3])
    emex = np.empty([len(selex),3])
    for k in range(3):
      mex[:,k] = f[selex,3+2*k]
      emex[:,k] = f[selex,3+2*k+1]
    for k in range(3):
      z = np.polyfit(attdB[sel],mm[sel,k][0],1)  # Linear fit
      yth = attdB*z[0]+z[1]
      #pl.subplot(2,1,1)
      fig = pl.figure(12)
      pl.errorbar(attdB,mm[:,k],yerr=emm[:,k],lw=2,label='Channel {0}'.format(k))
      pl.plot(attdB,yth,'y--')
      print 'Channel',k,', slope=',z[0],'V/dB'
      fig3 = pl.figure(13)
      sind = k+1
      sbp2 = pl.subplot(3,1,sind)
      sbp2.set_xscale('log')
      pl.errorbar(vin,mm[:,k],yerr=emm[:,k],lw=2,label='Channel {0}'.format(k))
      pl.errorbar(Vindaq,mex[:,k],yerr=emex[:,k],lw=2,label='ExtSin - Channel {0}'.format(k))
      pl.grid(True)
      if k==2:
    	pl.xlabel('Signal amplitude @ channel input [Vpp]')
      pl.ylabel('Mean output level [V]')
      pl.legend(loc='best')
  
  pl.figure(12)
  pl.title('Board {0} R{1}-{2}'.format(boardID,runstart,runstop))
  pl.xlabel('Quartz signal attenuation [dB]')
  pl.ylabel('Mean output level [V]')
  pl.legend(loc='best')
  pl.figure(12)
  pl.title('Board {0} R{1}-{2}'.format(boardID,runstart,runstop))
  
  
   
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
  
  if sys.argv[1] == "26":
    #runs=range(394,407)  # board 26 02/07/2018
    runs=range(470,496)  # board 12 02/07/2018
    
  if sys.argv[1] == "12":
    #runs=range(407,419)  # board 12 02/07/2018
    runs=range(444,470)  # board 12 02/07/2018
    
  if sys.argv[1] == "09":
    runs=range(3905,3918)  # board 09 27/08/2018

  if sys.argv[1] == "25":
    runs=range(4101,4114)  # board 25 27/08/2018

  if sys.argv[1] == "31":
    runs=range(4883,4895)  # board 25 27/08/2018


  #anaRuns(sys.argv[1],runs[0],runs[-1])
  anaRes(sys.argv[1],runs[0],runs[-1])

