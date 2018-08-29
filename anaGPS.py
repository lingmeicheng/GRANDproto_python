# -*- coding: utf-8 -*-
"""
Created on Fri Jun 08 10:19:36 2018

@author: Enya
"""

import numpy as np
import matplotlib.pyplot as plt
from pyproj import Proj, transform

ucenterfolder = "/home/martineau/GRAND/GRANDproto35/elec/tests/GPS/u-center/"
DISPLAY = 1

def readLog(filename):
  ubloxData = np.loadtxt(filename, dtype=str, delimiter='\n')
  lat, lon = [],[]
  ts = np.array([])
  for i in range(len(ubloxData)):
    line = ubloxData[i]
    myseg = line[1:6]   
    if myseg == "GPGLL": # Latitude / longitude and time values
        A = line.split(',')
        if A[2] == "N": # Interesting values lines
            a1 = A[1] #latitude 
            a2 = A[3] #longitude
            if len(a1)>1:
                b1 = float(a1)
            if len(a2)>1:
                b2 = float(a2)
            lat = np.append(lat, b1)
            lon = np.append(lon, b2)
            a3 = A[5] # time in hours/minutes/seconds
            t = (float(a3[0:2])*60*60+float(a3[2:4])*60+float(a3[4:6])) # Convert time in only seconds
            ts = np.append(ts, t)
  
  return lat,lon,ts            

def getXYpos(lat,lon):          
  inProj = Proj(init = 'epsg:3857') # Initial values coordinate system 
  Xt, Yt = [], []
  for i in range(len(lat)):
    x,y = inProj(lon[i]/100,lat[i]/100) # Convert lat and lon values in cartesians coordinates
    Xt = np.append(Xt, x)
    Yt = np.append(Yt, y)
  return Xt,Yt

# Read GPS log file for 1st antenna
[lat,lon,t1] = readLog(ucenterfolder+"COM3_180704_095810data.txt")
[Xt1,Yt1] = getXYpos(lat,lon)
Xtc1 = Xt1 - np.mean(Xt1) 
Ytc1 = Yt1 - np.mean(Yt1)
D1 = np.sqrt((Xtc1)**2 + (Ytc1)**2) # Distance between each measure and the [0;0] point
moy1 = np.mean(D1)
sigma1 = np.sqrt((np.sum((D1-moy1)**2))/len(D1)) # Standart deviation

# Read GPS log file for 2nd antenna
[lat,lon,t2] = readLog(ucenterfolder+"COM5_180704_095747data.txt")
[Xt2,Yt2] = getXYpos(lat,lon)
Xtc2 = Xt2 - np.mean(Xt2) 
Ytc2 = Yt2 - np.mean(Yt2)
D2 = np.sqrt((Xtc2)**2 + (Ytc2)**2) # Distance between each measure and the [0;0] point
moy2 = np.mean(D2)
sigma2 = np.sqrt((np.sum((D2-moy2)**2))/len(D2)) # Standart deviation

# Read time output data from anaTiming
dt = np.loadtxt('dPPS.txt')
dates_3 = np.loadtxt('dates.txt')


# Now build relative pos = f(t)
Dtot, Ttot = [], []
for i in range(len(t1)):
  j = np.argmin(abs(t2 - t1[i]))
  if abs(t2[j] - t1[i]) < 30:
      Ttot = np.append(Ttot, t2[j])
      Dtot = np.append(Dtot, np.sqrt((Xt1[i]-Xt2[j])**2 + (Yt1[i]-Yt2[j])**2))
Dtot2 = []
for i in range(len(dates_3)):
    h = np.argmin(abs(Ttot - dates_3[i]))
    Dtot2 = np.append(Dtot2, Dtot[h]) 

if DISPLAY:
  plt.scatter(dt, Dtot2)
  plt.xlabel('dt (s)')
  plt.ylabel('Dtot (m)')
  #plt.savefig('Dtot = f(dt).pdf', format='pdf')
  plt.show()

  plt.scatter(t1, D1)
  plt.title('dr1 = f(t1)')
  plt.xlabel('t1(s)')
  plt.ylabel('dr1 (m)')
  plt.savefig('dr1_t.pdf', format='pdf')
  plt.show()
  plt.scatter(t2, D2)
  plt.title('dr2 = f(t2)')
  plt.xlabel('t2(s)')
  plt.ylabel('dr2 (m)')
  plt.savefig('dr2_t.pdf', format='pdf')
  plt.show()

# Save resukts to files
np.savetxt('t1.txt', t1)
np.savetxt('t2.txt', t2)
np.savetxt('D1.txt', D1)
np.savetxt('D2.txt', D2)
np.savetxt('Dtot.txt', Dtot)


#raw_input()
#outProj = Proj(init = 'epsg:4326')
#x1,y1 = lon, lat
#x2, y2 = transform(inProj, outProj, x1, y1)
#print x2, y2      
#plt.scatter(x2, y2)
#plt.ylim(0.043576, 0.043577)
#plt.xlim(0.001987, 0.001989)
#plt.show() 
        

        
        
    
        
        
        
    
    
    
        




