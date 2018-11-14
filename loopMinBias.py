from anaMinBias import displayGalVar
from anaMinBias import loopRuns
import pylab as pl

ants = ["06","09","10","11","27","31"]
pl.ion()

for ant in ants:
  #loopRuns(ant,5968,6011)
  displayGalVar(ant)
  raw_input()
  pl.close('all')
