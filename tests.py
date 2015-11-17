from __future__ import division

import numpy as np

from fourier_tools import stft
from fourier_tools import plotstft

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

from pylab import plot
from pylab import subplot
from pylab import specgram
from pylab import show

import sys

from soundfiles_io import read_wav
from soundfiles_io import select_first_channel

from fourier_tools import stft
from fourier_tools import plotstft



# load sound
sound = "C:/repositories/ecoacoustics/CONAFOR-2015__1__20150304_093000.wav"

# you may import sounds in two ways, using scipy or the wave package
sample_rate_1, frames_1 = read_wav(sound,package="scipy")
#sample_rate_1, frames_1 = read_wav(sound,package="wave")

# sample rate of wave file
print (sample_rate_1)

# data from wave file is read as numpy array
print(np.shape(frames_1))

# plot spectrogram
plotstft(frames_1,sample_rate_1,numPoints=2048)

# get spectrogram data
fourier_data = stft(frames_1,numPoints=2048, frameSize=2**10)

# it is again a numpy array
print(np.shape(fourier_data))
