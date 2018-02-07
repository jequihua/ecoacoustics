# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 23:05:36 2018

@author: ttonaru
"""

import matplotlib.pyplot as plt
#import matplotlib.mlab as mlab
#
#from pylab import plot
#from pylab import subplot
#from pylab import specgram
#from pylab import show
#
#import sys
#
from soundfiles_io import read_wav
#from soundfiles_io import select_first_channel
#
#from fourier_tools import stft
#from fourier_tools import plotstft
# load sound
sound = ".\ccentralis_11.wav"
# read audiofile
fs, x = read_wav(sound,package="scipy")
# Detectio Threshold
theta = 0.4
# %%
# =============================================================================
# PYTHON IMPLEMENTATION
# =============================================================================
from vad_tools import PSDcompute,getvadparam,vadcompute,peakdetection,estnoisem
# computing of spectrogram as power spectrum density
Pxx,f,t,dt,df= PSDcompute(x,fs,'hamming',2048,0.75)
# adjust of parameters
vad,ne = getvadparam(dt,theta)
# noise estimation
Nxx1         = estnoisem(Pxx,dt,ne)
# implementation of a Voice Activity Detector 
Sx1,Sxx1       = vadcompute(Pxx,Nxx1,vad,ne)
# peak detection
det1 = peakdetection(Sx1,Sxx1,theta,t,f,vad,ne)
plt.plot(det1['DetectionSignal'][0])
#%%
# =============================================================================
# PYTHON + OCTAVE
# =============================================================================
#
#PYTHON + OCT2PY implementation
from oct2py import octave
octave.addpath('C:\Users\ttonaru\repos\ecoacoustics\batcall_detection')
# computing of spectrogram as power spectrum density
Pxx,f,t,dt,df= PSDcompute(x,fs,'hamming',2048,0.75)
# adjust of parameters
prm = octave.getvadparam(dt,theta)
vad = prm['vad']
ne  = prm['ne']
# noise estimation
Nxx2 = octave.estnoisems(Pxx,dt,ne)
# implementation of a Voice Activity Detector 
LR  = octave.Lratio(Pxx,Nxx2,vad,ne)
Sx2  = LR['Sx']
Sxx2 = LR['Sxx']
# peak detection
det2 = octave.signaldetect(Sx2,Sxx2,theta,t,f,prm)
plt.plot(det2['DetectionSignal'][0])



