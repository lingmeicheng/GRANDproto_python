import os
import time
import sys
import math

import numpy as np
import pylab as pl

print sys.argv
boardid = sys.argv[1]
folder = sys.argv[2]

fname = 'freq_b' + boardid + '.txt' 
fname = folder+'/'+fname
count = -1
f = open(fname)
freq = range(20,115,5)  # taken from testfreq.sh

m = np.zeros([np.size(freq),3])
sd = np.zeros([np.size(freq),3])
for line in iter(f):
  subline = line.split('= ')
  if np.size(subline)==3:
    ch  = int(subline[0].split(' ')[1])
    m[count,ch] = float(subline[1].split(' ')[0])
    sd[count,ch] = float(subline[2].split(' ')[0])

  elif np.size(subline)==1:
    runid = subline[0].split(' ')[6]
    count = count+1
f.close()

chl = ['X','Y','Z']
c = ['b','r','g']
pl.figure(1)
for k in range(3):
  pl.errorbar(freq,m[:,k],yerr=sd[:,k],color=c[k],linewidth=2,label=chl[k])
  
pl.grid(True)
pl.xlabel('Frequency (MHz)')
pl.ylabel('Output (V)')
pl.title('Freq calib board {0} - Input = Sin 100mVpp-30dB-splitter'.format(boardid))  
pl.legend(loc='best')
#pl.show()  

pl.savefig(folder+'/freq.png')
