
import wave

from scipy.io.wavfile import read
from scipy.io.wavfile import write

import numpy as np


def read_wav(file_name,package="scipy"):
	'''
	file_name: file_name or path + file_name

	package: which python package to use for reading the file
	
	'''

	if (package=="scipy"):
		fs, frames = read(file_name)
	elif (package=="wave"):
		read_wav = wave.open(file_name,'r')
		fs = read_wav.getframerate()
		read_wav = read_wav.readframes(read_wav.getnframes())
		frames = np.fromstring(read_wav, 'Int16')

	return (fs,frames)

def select_first_channel(frames,package="scipy"):

	'''
	frames: output produced by read_wav
	
	'''
	if (package=="scipy"):
		signal = frames[:,0]
	elif (package=="wave"):
		signal = frames[0::2]

	return (signal)





