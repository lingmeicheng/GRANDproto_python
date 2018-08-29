import numpy as np
from anaMinBias import loopEvents

runs = np.arange(1200,1400)

for run in runs:
  print 'Launching run',run
  loopEvents(str(run),"05")
  
  
  
