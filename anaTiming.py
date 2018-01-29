import os
import time
import sys
import mx.DateTime
import math

import numpy as np
import pylab as pl

datadir = "/home/martineau/GRAND/GRANDproto35/data/tests/"

def loopEvents(filename):

   datafile = datadir+filename
   print 'Scanning',datafile

   with open(datafile,"r") as f:
   	   evts = f.read().split('-----------------')

   print 'Number of events:',len(evts)-1
   time.sleep(1)
   date = []
   board = np.zeros(shape=(np.size(evts)))
   TS2 = np.zeros(shape=(np.size(evts)))
   TS1PPS = np.zeros(shape=(np.size(evts)))
   TS1Trig = np.zeros(shape=(np.size(evts)))
   SSS = np.zeros(shape=(np.size(evts)))
   EvtId = np.zeros(shape=(np.size(evts)))
   TrigPattern = np.zeros(shape=(np.size(evts)))
   trigtime = np.zeros(shape=(np.size(evts)))
   data = list()
   j = 0

   for i in range(1001,np.size(evts)):  # Skip at least first elemt (empty)
   #for i in range(1001,2001):  # Skip first elemt which is empty
     if float(i)/100 == int(i/100):
     	  print 'Event ',i

 
     evt = evts[i]
     evtsplit = evt.split('\n')
     if np.size(evtsplit)>8:   # Event is of normal size
       #date.append(evtsplit[1])
       d = mx.DateTime.DateTimeFrom(evtsplit[1])  # Better parser than standatetime module
       date.append(d)
       IP = evtsplit[2][3:]
       if IP=='192.168.1.101':
         board[j] = 01
       if IP=='192.168.1.102':
         board[j] = 02

       TS2[j]=int(evtsplit[3][4:])  # time elapsed since last PPS (125MHz clock <=> ounter)
       tt=int(evtsplit[4][11:])  # phase in 8ns slot fr trigger
       TS1Trig[j] = get_1stone(hex(tt))
       tpps=int(evtsplit[5][7:])
       TS1PPS[j]=get_1stone(hex(tpps))  # phase in 8ns slot for PPS
       SSS[j]=int(evtsplit[6][4:])  # Elapsed seconds since start
       EvtId[j] = int(evtsplit[7][3:])
       TrigPattern[j] = int(evtsplit[8][12:])
       # Data
       raw=evtsplit[9:][:]  #raw data
       raw2 = raw[0].split(" ") # Cut raw data list into samples
       raw2 = raw2[0:np.size(raw2)-1]	# Remove last element (empty)
       hraw2 = [hex(int(a)) for a in raw2]  # TRansfer back to hexadecimal
       draw = [twos_comp(int(a,16), 12) for a in hraw2] #2s complements

       nsamples = len(draw)/4  # Separate data to each channel
       offset = nsamples/2.0
       #print nsamples,"samples per channel --> offset = ",offset
       data.append(np.reshape(draw,(4,nsamples))) # Write to data list
     else:
       print 'Error! Empty event',i
       date.append("")
       data.append("")

     j = j+1

   # Build trig time
   date = np.array(date)
   boards = set(board[np.where(board>0)])
   print 'Boards in run:',list(boards)
   print 'Run start:', date[0]
   
   sel1 = np.where(board == list(boards)[0])  # Use same PPS second for all cards
   for id in boards:
     sel = np.where(board == id)
     date_end = date[sel[0][-1]-1]
     print 'Run stop:',date_end,'for board',id,' (',np.size(sel),'measurements)'
     #pl.figure(12)
     #pl.hist(TS2[sel])
     #pl.xlabel("Distrib of TS2 (8ns step counter)")


     # Grab MaxCoarse info from SLC run.
     if id == 1:
       [d,maxCoarseRaw] = getMaxCoarse("S22_b01.data")
     if id == 2:
       [d,maxCoarseRaw] = getMaxCoarse("S22_b02.data")
     d = np.array(d)
     maxCoarse = np.zeros(np.shape(sel))
     for k in range(np.size(sel)):
     	idift = np.argmin(np.abs(date[k]-d))
	maxCoarse[0,k] = maxCoarseRaw[idift]
	
     # Now compute trig time
     #cor=1
     #cor=125e6/float(max(TS2[sel]))
     cor=125e6/(maxCoarse+1)
     print 'Correction factor for 125MHz clock for board',id,':',cor
     trigtime[sel] = SSS[sel1] +(TS2[sel]*4+TS1PPS[sel]-TS1Trig[sel])*2e-9*cor  #s. Use same SSS for both cards. Warning: requires 100% trigger + odd nb of events
     
     if id == 1:
       newS1 = np.where(np.diff(SSS[sel])>0)[[0][-1]]+1  # index of 1st event in new second for board 1
     if id == 2:
       newS2 = np.where(np.diff(SSS[sel])>0)[[0][-1]]+1  # index of 1st event in new second for board 2

   trigtime1 = trigtime[np.where(board == list(boards)[0])]
   trigtime2 = trigtime[np.where(board == list(boards)[1])]
   #TS21 = TS2[np.where(board == list(boards)[0])]  # This is the clock count for the 1st event in the new second
   tdiff = (trigtime2-trigtime1)*1e9  # Delta t in ns
   dPPS = (trigtime2[newS2]-trigtime1[newS1])*1e9  # Delta t @ new PPS in ns

   tdiff_b1 = np.zeros(shape=(np.size(trigtime1)))
   for i in range(len(trigtime1)-1):
   	tdiff_b1[i]=(trigtime1[i+1]-trigtime1[i])*1e3
   tdiff_b2 = np.zeros(shape=(np.size(trigtime2)))
   for i in range(len(trigtime2)-1):
   	tdiff_b2[i]=(trigtime2[i+1]-trigtime2[i])*1e3

   tdiff_b1 = tdiff_b1[tdiff_b1>0]
   tdiff_b2 = tdiff_b2[tdiff_b2>0]


   print 'Time delay beween consecutive events for board 1:',np.mean(tdiff_b1),'ms'
   print 'Time delay beween consecutive events for board 2:',np.mean(tdiff_b2),'ms'

   sel = np.where(np.abs(tdiff)<300)[[0][-1]]
   print 'Delta Trig Time = ',np.mean(tdiff[sel]),'+-',np.std(tdiff[sel]),'ns'
   print 'Delta Trig Time @ new PPS = ',np.mean(dPPS),'+-',np.std(dPPS),'ns'
   pl.figure(1)
   pl.plot(trigtime1[1:],label='Board 1')
   pl.hold(True)
   pl.grid(True)
   pl.plot(trigtime2[1:],'r',label='Board 2')
   pl.xlabel('Event nb')
   pl.ylabel('Trig time [s]')
   pl.legend()

   pl.figure(21)
   pl.subplot(311)
   pl.plot(trigtime1,tdiff)
   pl.plot(trigtime1[newS1],dPPS,marker="o")
   pl.xlabel('Time [s]')
   pl.ylabel('trigtime_b2-trigtime_b1 [ns]')
   pl.grid(True)
   pl.xlim(min(trigtime1),max(trigtime1))
   pl.subplot(312)
   pl.hist(tdiff[sel],100)
   pl.xlabel('trigtime_b2-trigtime_b1 [ns]')
   pl.subplot(313)
   pl.hist(dPPS,100)
   pl.xlabel('trigtime_b2-trigtime_b1 @ newPPS [ns]')

   pl.figure(18)
   pl.subplot(221)
   pl.plot(tdiff_b1[1:])
   pl.xlabel('Event nb')
   pl.ylabel('Consecutive Delta Trigtime [ms]')
   pl.title("Board b1")
   pl.subplot(222)
   pl.title("Board b1")
   pl.hist(tdiff_b1[1:],1000)
   pl.xlabel('Consecutive Delta Trigtime [ms]')
   pl.subplot(223)
   pl.title("Board b2")
   pl.plot(tdiff_b2[1:])
   pl.xlabel('Event nb')
   pl.ylabel('Consecutive Delta Trigtime [ms]')
   pl.subplot(224)
   pl.title("Board b2")
   pl.hist(tdiff_b2[1:],1000)
   pl.xlabel('Consecutive Delta Trigtime [ms]')
   pl.show()


def getMaxCoarse(filename):
# retrieve maxCoarse info from bslow control data
    datafile = datadir+filename
    print 'Scanning',datafile

    with open(datafile,"r") as f:
    	   evts = f.read().split('-----------------')

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
