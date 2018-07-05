import numpy as np
from anaMinBias import loopEvents

runs = np.arange(1200,1400)

for run in runs:
  print run
  loopEvents(str(run),"05")
  
  
  
