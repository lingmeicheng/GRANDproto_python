import os
import time
import sys
import math

import numpy as np
import pylab as pl

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 18}

pl.rc('font', **font)

consigne = np.array([100, 1000, 5000, 6000, 8000, 10000, 12000, 14000, 16000, 18000, 20000])  #Hz
trate = 2*consigne*1e-3
dur = np.array([7.58, 12.89, 2.77, 2.8, 2.41, 2.7, 2.02, 2.27, 2.2, 1.96, 1.64])  # s
nev = np.array([1517, 25772, 27780, 33687, 38635, 53990, 41358, 46083, 48231, 38024, 36167])
rate = nev/dur*1e-3
expev = 2*consigne*dur
eff = nev/expev

cons45 = np.array([6000, 12000, 14000, 16000, 18000, 20000])
trate45 = 2*cons45*1e-3
dur45 = np.array([1.13, 2.17, 1.97, 2.00, 1.1, 1.78,])
nev45 = np.array([13661, 52190, 53795, 63604, 39432, 53195])
rate45 = nev45/dur45*1e-3
expev = 2*cons45*dur45
eff45 = nev45/expev


pl.figure(1)
pl.plot(trate,eff,'s',ms=8)
pl.plot(trate45,eff45,'o',ms=8)
pl.grid(True)
pl.xlabel('Trig rate [kHz]')
pl.ylabel('Recorded fraction [%]')
axes = pl.gca()
axes.set_xlim([0,42])

pl.figure(2)
pl.plot(trate,rate,'s',ms=8)
pl.plot(trate45,rate45,'o',ms=8)
pl.grid(True)
pl.xlabel('Trig rate [kHz]')
pl.ylabel('Record rate [kHz]')
axes = pl.gca()
axes.set_xlim([0,42])

