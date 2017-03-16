import os
import time
import sys
import math

import numpy as np
import pylab as pl

def loop():
   datafile = '../data/slc_remote.txt'
   print 'Scanning',datafile

   with open(datafile,"r") as f:
   	   evts = f.read().split('-----------------')

   nevts=len(evts)-1
   print 'Number of events:',nevts
   
   date = []
   board = np.zeros(shape=(np.size(evts)),dtype = np.int32)
   VPower = np.zeros(shape=(nevts,6))
   Th = np.zeros(shape=(nevts,3,2))
   Temp = np.zeros(shape=(np.size(evts)))
   TrigRate = np.zeros(shape=(nevts,7))
   j = 0
   for i in range(j,nevts):
   	evt = evts[i]
	print evt
	print len(evt)
   	if len(evt)>300:
	     print 'SLC word'
	     
	     evtsplit = evt.split('\n')
	     date.append(evtsplit[1])
	     IP = evtsplit[2][3:]
	     board[j] = int(IP[-2:]);
	     for k in range(0,5):
   	       VPower[j,k]=int(evtsplit[k+3][8:])  
	     print VPower[j,:]
	     for k in range(0,2):
	       thsplit = evt.split(':')
	       Th[j,k,1] = int(evtsplit[k+9][12:1]) 
	     j = j+1

      	else:
   	     print 'ACK word, skip it.'
   	     j = i
 
if __name__ == '__main__':
     loop()
