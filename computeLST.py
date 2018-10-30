import time
import ephem
import datetime
import sys

# Gansu 
longi = "94.10"
lat = "39.35"

y = int(sys.argv[1])
m = int(sys.argv[2])
d = int(sys.argv[3])
h, mn, s = 0, 0 , 0
if len(sys.argv)>4:
  h = int(sys.argv[4])
  if len(sys.argv)>5:
    mn = int(sys.argv[5])
    if len(sys.argv)>6:
      s = int(sys.argv[6])
  
site = ephem.Observer();
site.date = datetime.datetime(y,m,d,h,mn,s)   # GMT time
site.long = ephem.degrees(longi)
site.lat = ephem.degrees(lat)
site.elevation = 1500;

print "Site: ({0}N,{1}E)".format(longi,lat)
print "GMT time: {0}/{1}/{2} {3}:{4}:{5}".format(y,m,d,h,mn,s)
print "Sideral time:",site.sidereal_time()
