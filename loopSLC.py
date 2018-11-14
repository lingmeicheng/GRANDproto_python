from anaSLC import displaySLC
from anaSLC import loopSLCEvents
import pylab as pl

ants = ["06","09","10","11","27","31"]
pl.ion()

for ant in ants:
  loopSLCEvents(ant,5968,6011)
  #displaySLC(ant)
  #raw_input()
  #pl.close('all')
