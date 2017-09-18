#!/home/sandra/anaconda3/bin/ipython

import os
CC=0

AnalysisType=1
ALL3DANTS=[1, 2, 3, 4, 5, 6]
ALL3DSCINTS=[111, 112, 113, 114, 115, 116]
ALL3DDETS=ALL3DANTS+ALL3DSCINTS

RAW_DATA_PATH='/home/sandra/data/'
if CC:
    tmpdir = os.getenv("TMPDIR")
    RAW_DATA_PATH=tmpdir+'/'


STD_PATH='/home/sandra/Python/GRANDprotoSoft/'
if CC:
    STD_PATH='/afs/in2p3.fr/throng/trend/soft/rec/pythonGRANDproto/'

PLOT_PATH=STD_PATH+'plot/'

if AnalysisType==0:
    DST_DATA_PATH=STD_PATH+'dstant/'
if AnalysisType==1:
    DST_DATA_PATH=STD_PATH+'dsthyb/'
if AnalysisType==2:
    DST_DATA_PATH=STD_PATH+'dstsci/'


#DST_DATA_PATH=STD_PATH+'dstant/'
TEXT_PATH = STD_PATH+'txt/'
CXX_PATH='/home/sandra/TREND/cxx/'
LOG_PATH='/home/sandra/monitoringDAQ/'
ibuffs=1024
ibuff=1024
ibufft=4
TSAMPLING=5e-9 #sec
FSAMPLING=200e6 #Hz
TriggerTimeSpan=60. #s
C0=2.99792458e8 #m/s
CORREL=0
ErrorTrig=2
REFWE=0
REFSN=0
REFALT=2632.
#SCALE=0.82/74 #trend ou gproto?
SCALE=280e-3/39
FREQMIN=30e6 #trend ou gproto?
FREQMAX=100e6
ConsCoincCompleteRejection=1
threshold=5.0
granularity_box=40e-9 #sec
granularity_sampling=15e-9 #sec
max_out=0
max_out_ToT=300e-9 #sec
max_block_ToT=400e-9 #sec
noise_mode=1
QuickReject=1
RawFilterTag=1
max_pretrig_ToT=200e-9
