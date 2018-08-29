import os
import time
import sys
import mx.DateTime
import math

import numpy as np
import pylab as pl

# Local module. Do ./setup.py build & ./setup.py install if not loadable [Junhua Sept 27 2017]
import pyef 

#datadir = "/home-local/GRANDproto/GRANDproto/data"
datadir = "/home/martineau/GRAND/GRANDproto35/elec/tests/GPS/data"

def loopEvents(RUNID,ID1,ID2):
   filename = 'R'+str(RUNID)+'_b'+str(ID2)+'.data.txt'
   
   datafile = datadir+'/'+filename
   print 'Scanning',datafile

   with open(datafile,"r") as f:
   	   evts = f.read().split('-----------------')

   
   #########################################
   ### Set nb of events
   #########################################
   minN = 1
   #maxN = 1000
   maxN = np.size(evts)
   if (maxN-minN)%2 == 1: 
     print "Warning... Odd number of events. Cutting last one to build pairs."
     maxN = maxN-1
   dimarray = maxN-minN+1

   print 'Number of events:',len(evts)-1
   time.sleep(1)
   date = []
   dates = []
   board = np.zeros(shape=(dimarray))
   TS2 = np.zeros(shape=(dimarray))
   TS1PPS = np.zeros(shape=(dimarray))
   TS1Trig = np.zeros(shape=(dimarray))
   SSS = np.zeros(shape=(dimarray))
   EvtId = np.zeros(shape=(dimarray))
   TrigPattern = np.zeros(shape=(dimarray))
   trigtime = np.zeros(shape=(dimarray))
   data = list()
   j = 0
         
   for i in range(minN,maxN):  # Here toggle+0/-1 if crashes at trigtime[sel] = SSS[sel1] +(TS2[sel]*4+TS1PPS[sel]-TS1Trig[sel])*2e-9*cor (ugly fix) 
     
     if float(i)/100 == int(i/100):
     	  print 'Event ',i

 
     evt = evts[i]
     evtsplit = evt.split('\n')
     if np.size(evtsplit)>8:   # Event is of normal size
       #date.append(evtsplit[1])
       d = mx.DateTime.DateTimeFrom(evtsplit[1])  # Better parser than standatetime module
       evt1 = evtsplit[1].split(' ')
       timefield = evt1[-3]
       ev1 = (float(timefield[0:2])*60*60+float(timefield[3:5])*60+float(timefield[6:]))
       dates = np.append(dates, ev1)
       
       date.append(d)
       IP=int(evtsplit[2].split('.')[-1])
       if IP==100+int(ID1):
         board[j] = 01
       if IP==100+int(ID2):
         board[j] = 02
       if j>0 and board[j]==board[j-1]:
	  print i,j,board[j-1],board[j],int(evtsplit[6][4:])
	  print "## Error! 2 consecutive events from same board! Erase 1st one."
          date.pop()
	  j = j-1
	  #raw_input()
	  #continue
	    
       TS2[j]=int(evtsplit[3][4:])  # time elapsed since last PPS (125MHz clock <=> ounter)
       #print j%2,": board",int(board[j]),"IP=",IP,"TS2=",TS2[j]
       if j%2==1 and abs(TS2[j]-TS2[j-1])>1000:
          print "## Error! Mismatch in TS2 timings...",TS2[j-1],TS2[j],". Deleting this pair."
          j = j-1
          continue
	  
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
   
   if board[-1]==board[0]:
     print "Skipping last (unpaired) event"
     board = board[:-1]
     
   # Build trig time
   date = np.array(date)
   boards = set(board[np.where(board>0)])
   print 'Boards in run:',list(boards)
   print 'Run start:', date[0]
   
   sel1 = np.where(board == list(boards)[0])  # Compute time with same PPS second for all cards
   for id in boards:
     sel = np.where(board == id)
     date_end = date[sel[0][-1]-1]
     print 'Run stop:',date_end,'for board',id,' (',np.size(sel),'measurements)'
     #pl.figure(12)
     #pl.hist(TS2[sel])
     #pl.xlabel("Distrib of TS2 (8ns step counter)")


     # Grab MaxCoarse info from SLC run.
     if id == 1:
       [d,maxCoarseRaw] = getMaxCoarse('S'+str(RUNID)+'_b'+str(ID1)+'.data.txt')
       #[d,maxCoarseRaw] = getMaxCoarse('S1110_b'+str(ID1)+'.data.txt')
     if id == 2:
       [d,maxCoarseRaw] = getMaxCoarse('S'+str(RUNID)+'_b'+str(ID2)+'.data.txt')
       #[d,maxCoarseRaw] = getMaxCoarse('S1110_b'+str(ID2)+'.data.txt')
     
     d = np.array(d)
     if np.size(d)==1:
        print "No maxCoarse data ==> no correction."
        cor = 1
     else:
       maxCoarse = np.zeros((np.size(sel),1))
       for k in range(np.size(sel)): 
     	 idift = np.argmin(np.abs(date[k]-d))
	 maxCoarse[k] = maxCoarseRaw[idift]
       cor=125e6/(maxCoarse+1)
       pl.figure(13)
       if id==1:
         pl.plot(maxCoarse,'+b',label='board 1')
       if id==2:
	 pl.plot(maxCoarse,'+r',label='board 2')
       pl.xlabel('Event nb')
       pl.ylabel('MaxCoarse')
       pl.legend()
       #print 'Correction factor for 125MHz clock for board',id,':',cor
     
     trigtime[sel] = SSS[sel1] +(TS2[sel]*4+TS1PPS[sel]-TS1Trig[sel])*2e-9*cor  #s. Use same SSS for both cards. Warning: requires 100% trigger + odd nb of events
     if id == 1:
       newS1 = np.where(np.diff(SSS[sel])>0)[[0][-1]]+1  # index of 1st event in new second for board 1
     
     #if id == 2:
     #  newS2 = np.where(np.diff(SSS[sel])>0)[[0][-1]]+1  # index of 1st event in new second for board 2
     #  maxIndex = np.min([np.size(newS1),np.size(newS2)])
     #  print maxIndex
     #  newS1 = newS1[0:maxIndex]
     #  newS2 = newS2[0:maxIndex]
     #  print newS1
     #  print newS2

   dates_2 = dates[newS1]
   trigtime1 = trigtime[np.where(board == list(boards)[0])]
   trigtime2 = trigtime[np.where(board == list(boards)[1])]
   #TS21 = TS2[np.where(board == list(boards)[0])]  # This is the clock count for the 1st event in the new second
   tdiff = (trigtime2-trigtime1)*1e9  # Delta t in ns
   dPPS = (trigtime2[newS1]-trigtime1[newS1])*1e9  # Delta t @ new PPS in ns

   tdiff_b1 = np.zeros(shape=(np.size(trigtime1)))
   for i in range(len(trigtime1)-1):
   	tdiff_b1[i]=(trigtime1[i+1]-trigtime1[i])*1e3
   tdiff_b2 = np.zeros(shape=(np.size(trigtime2)))
   for i in range(len(trigtime2)-1):
   	tdiff_b2[i]=(trigtime2[i+1]-trigtime2[i])*1e3

   #tdiff_b1 = tdiff_b1[tdiff_b1>0]
   #tdiff_b2 = tdiff_b2[tdiff_b2>0]

   np.savetxt('dates.txt', dates_2)
   np.savetxt('dPPS.txt', dPPS)

   selb1 = np.intersect1d(np.where(tdiff_b1<np.mean(tdiff_b1)+np.std(tdiff_b1)),np.where(tdiff_b1>np.mean(tdiff_b1)-np.std(tdiff_b1)))
   selb2 = np.intersect1d(np.where(tdiff_b2<np.mean(tdiff_b2)+np.std(tdiff_b2)),np.where(tdiff_b2>np.mean(tdiff_b2)-np.std(tdiff_b2)))
   dif1 = np.mean(tdiff_b1[selb1])
   dif2 = np.mean(tdiff_b2[selb2])
   #print 'Time delay beween consecutive events for board 1:',dif1,'ms'
   #print 'Time delay beween consecutive events for board 2:',dif2,'ms'

   sel = np.where(np.abs(tdiff)<200)[[0][-1]]
   selpps = np.where(np.abs(dPPS)<200)[[0][-1]]
   print 'Delta Trig Time = ',np.mean(tdiff[sel]),'+-',np.std(tdiff[sel]),'ns'
   print 'Delta Trig Time @ new PPS = ',np.mean(dPPS[selpps]),'+-',np.std(dPPS[selpps]),'ns'
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
   pl.hist(dPPS[selpps],100)
   pl.xlabel('trigtime_b2-trigtime_b1 @ newPPS [ns]')
   
   selb1 = np.intersect1d(np.where(tdiff_b1<dif1+0.001),np.where(tdiff_b1>dif1-0.001))
   selb2 = np.intersect1d(np.where(tdiff_b2<dif2+0.001),np.where(tdiff_b2>dif2-0.001))
   pl.figure(18)
   pl.subplot(221)
   pl.plot(tdiff_b1[selb1]*1e6)
   pl.xlabel('Event nb')
   pl.ylabel('Consecutive Delta Trigtime (ns)')
   pl.title("Board b1")
   pl.subplot(222)
   pl.title("Board b1")
   pl.hist(tdiff_b1[selb1]*1e6,1000)
   pl.xlabel('Consecutive Delta Trigtime (ns)')
   pl.subplot(223)
   pl.title("Board b2")
   pl.plot(tdiff_b2[selb2]*1e6)
   pl.xlabel('Event nb')
   pl.ylabel('Consecutive Delta Trigtime (ns)')
   pl.subplot(224)
   pl.title("Board b2")
   pl.hist(tdiff_b2[selb2]*1e6,1000)
   pl.xlabel('Consecutive Delta Trigtime (ns)')
   print 'Time delay beween consecutive events for board 1:',np.mean(tdiff_b1[selb1]),'ms +-',np.std(tdiff_b1[selb1]*1e6),'ns'
   print 'Time delay beween consecutive events for board 2:',np.mean(tdiff_b2[selb2]),'ms +-',np.std(tdiff_b2[selb2]*1e6),'ns'
   pl.show()


def getMaxCoarse(filename):
# retrieve maxCoarse info from bslow control data
    datafile = datadir+'/'+filename
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
     loopEvents(sys.argv[1],sys.argv[2],sys.argv[3])
