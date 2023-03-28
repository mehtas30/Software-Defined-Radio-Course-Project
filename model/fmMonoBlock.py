# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada

import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy import signal
import numpy as np
import math
import sys

from fmSupportLib import fmDemodArctan, fmPlotPSD

def lp_impulse_response_coeff(Fc, Fs, N_taps):
	normCutoff = Fc/(Fs/2)

	coefficients = np.zeros(N_taps)

	for i in range(N_taps):
		if (i == ((N_taps-1)/2)):
			coefficients[i] = normCutoff

		else:
			# sinc = sin(x)/x
			denominator = math.pi * normCutoff * (i - ((N_taps-1) / 2))
			numerator = math.sin(denominator)
			

			# normalize coefficients
			coefficients[i] = normCutoff * (numerator/denominator)

		# window the coefficients
		coefficients[i] = coefficients[i] * ((math.sin((i*math.pi) / N_taps)) ** 2)

	return coefficients

def lp_filter(coefficients, data, state):

	state_len = len(state)
	data_len = len(data)
	coeff_len = len(coefficients)
	filtered_data = np.zeros(data_len)
	
	data = np.concatenate([state,data])

	# discrete convolution
	for n in range(data_len):
		for k in range(len(coefficients)):
			filtered_data[n] += coefficients[coeff_len - k -1] * data[n+k]

	# current unfiltered block is next block's filter state
	filter_state = data[-state_len:]
	
	return filtered_data, filter_state

def myDemod(i_ds, q_ds, p_i=0, p_q=0):
	prevI, prevQ = p_i, p_q
	demod = np.array([])

	for i in range (len(i_ds)):
		currentI = i_ds[i]
		currentQ = q_ds[i]

		derivI = currentI - prevI
		derivQ = currentQ - prevQ
		
		denominator = currentI**2 + currentQ**2

		if (denominator == 0):
			demod = np.append(demod, 0)
		else:
			numerator = currentI*derivQ - currentQ*derivI
			demod = np.append(demod, numerator/denominator)

		prevI = currentI
		prevQ = currentQ

	return demod, prevI, prevQ
	
def upsample(data, up_factor):
	if up_factor == 1:
		return data
		
	resampled_data = np.zeros(len(data) * up_factor)
	for i in range(len(data)):
		resampled_data[i * up_factor] = data[i]
	return resampled_data
	
def downsample(data, down_factor):
	return data[::down_factor]
			
	
	
mode = 0

rf_Fs = 2.4e6
rf_Fc = 100e3
rf_taps = 151
rf_decim = 10

if_Fs = 240e3

audio_Fs = 48e3
audio_Fc = 16e3
audio_taps = 101
audio_decim = 5
audio_interp = 1

'''
mode 0:
2.4MS/s ---10 down---> 240kS/s ---5 down---> 48kS/s

mode 1:
1.152MS/s ---4 down---> 288kS/s ---4 down---> 48kS/s

mode 2:
2.4MS/s ---10 down---> 240kS/s ---147 up, 800 down---> 44.1kS/s

mode 3:
2.304MS/s ---9 down---> 256kS/s ---441 up, 2560 down---> 44.1kS/s
'''

if __name__ == "__main__":
	
	if mode == 0:
		rf_Fs = 2.4e6
		rf_decim = 10
		if_Fs = 240e3
		audio_Fs = 48e3
		audio_decim = 5
	elif mode == 1:
		rf_Fs = 1.152e6
		rf_decim = 4
		if_Fs = 288e3
		audio_Fs = 48e3
		audio_decim = 4
	elif mode == 2:
		rf_Fs = 2.4e6
		rf_decim = 10
		if_Fs = 240e3
		audio_Fs = 44.1e3
		audio_decim = 800
		audio_interp = 147
		
		audio_taps *= audio_interp
		if_Fs *= audio_interp
	elif mode == 3:
		rf_Fs = 2.304e6
		rf_decim = 9
		if_Fs = 256e3
		audio_Fs = 44.1e3
		audio_decim = 2560
		audio_interp = 441

		audio_taps *= audio_interp
		if_Fs *= audio_interp
	else:
		print(f'Invalid operating mode!')
		sys.exit()

	# read the raw IQ data from the recorded file
	# IQ data is assumed to be in 8-bits unsigned (and interleaved)
	in_fname = "../data/samples_u8.raw"
	raw_data = np.fromfile(in_fname, dtype='uint8')
	print("Read raw RF data from \"" + in_fname + "\" in unsigned 8-bit format")
	# IQ data is normalized between -1 and +1 in 32-bit float format
	iq_data = (np.float32(raw_data) - 128.0)/128.0
	print("Reformatted raw RF data to 32-bit float format (" + str(iq_data.size * iq_data.itemsize) + " bytes)")

	# set up the subfigures for plotting
	# subfig_height = np.array([0.8, 2, 1.6]) # relative heights of the subfigures
	# plt.rc('figure', figsize=(7.5, 7.5))	# the size of the entire figure
	# fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, gridspec_kw={'height_ratios': subfig_height})
	# fig.subplots_adjust(hspace = .6)
	
	print(f'Operating in mode {mode}\n\n\
			rf_Fs = {rf_Fs}\n\
			rf_Fc = {rf_Fc}\n\
			rf_taps = {rf_taps}\n\
			rf_decim = {rf_decim}\n\n\
			if_Fs = {if_Fs}\n\n\
			audio_Fs = {audio_Fs}\n\
			audio_Fc = {audio_Fc}\n\
			audio_taps = {audio_taps}\n\
			audio_decim = {audio_decim}\n\
			audio_interp = {audio_interp}\n')

	# select a block_size that is a multiple of KB
	# and a multiple of decimation factors
	block_size = 128 * rf_decim * audio_decim* 2
	block_count = 0

	# coefficients for IQ -> IF LPFs, Fc = 100kHz
	rf_coeff = lp_impulse_response_coeff(rf_Fc, rf_Fs, rf_taps)

	# coefficients for IF -> audio LPF, Fc = 16kHz
	audio_coeff = lp_impulse_response_coeff(audio_Fc, if_Fs, audio_taps)
	audio_coeff *= audio_interp

	# state-saving for extracting FM band (Fc = 100kHz)
	state_i_lpf_100k = np.zeros(rf_taps-1)
	state_q_lpf_100k = np.zeros(rf_taps-1)
	# state-saving for FM demodulation
	prevI = 0
	prevQ = 0
	# state-saving for extracting mono audio (Fc = 16kHz)
	audio_state = np.zeros(audio_taps-1)

	# audio buffer that stores all the audio blocks
	audio_data = np.array([]) # used to concatenate filtered blocks (audio data)

	# if the number of samples in the last block is less than the block size
	# it is fine to ignore the last few samples from the raw IQ file
	while (block_count+1)*block_size < len(iq_data):

		# if you wish to have shorter runtimes while troubleshooting
		# you can control the above loop exit condition as you see fit
		print(f'Processing block {block_count}')

		# filter to extract the FM channel (I samples are even, Q samples are odd)
		i_filt, state_i_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
				iq_data[(block_count)*block_size:(block_count+1)*block_size:2],
				zi=state_i_lpf_100k)
		q_filt, state_q_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
				iq_data[(block_count)*block_size+1:(block_count+1)*block_size:2],
				zi=state_q_lpf_100k)

		# downsample the I/Q data from the FM channel
		i_ds = downsample(i_filt, rf_decim)
		q_ds = downsample(q_filt, rf_decim)
											
		
		# FM demodulator
		# you will need to implement your own FM demodulation based on:
		# https://www.embedded.com/dsp-tricks-frequency-demodulation-algorithms/
		# see more comments on fmSupportLib.py - take particular notice that
		# you MUST have also "custom" state-saving for your own FM demodulator
		fm_demod, prevI, prevQ = myDemod(i_ds, q_ds, prevI, prevQ)
		
		# if (block_count == 0):
			# print('fm_demod=')
			# for i in range(30):
				# print(fm_demod[i])
		
		# extract mono audio data through filtering, and downsample to 48kS/s
		fm_demod_us = upsample(fm_demod, audio_interp)

		# extract the mono audio data through filtering
		audio_filt, audio_state = signal.lfilter(audio_coeff, 1.0, fm_demod_us, zi=audio_state)

		# downsample audio data
		audio_block = downsample(audio_filt, audio_decim)
		
		if (block_count == 5):
			print('audio_block=')
			for i in range(30):
				print(audio_block[i])

		# concatenate the most recently processed audio_block
		# to the previous blocks stored already in audio_data
		audio_data = np.concatenate((audio_data, audio_block))

		# to save runtime select the range of blocks to log data
		# this includes both saving binary files as well plotting PSD
		# below we assume we want to plot for graphs for blocks 10 and 11
		#if block_count == 0:

		 	# plot PSD of selected block after FM demodulation
		 	#ax0.clear()

		 	#fmPlotPSD(ax0, fm_demod, (rf_Fs/rf_decim)/1e3, subfig_height[0], \
		 	#		'Demodulated FM (block ' + str(block_count) + ')')
		 	# output binary file name (where samples are written from Python)
		 	#fm_demod_fname = "../data/fm_demod_" + str(block_count) + ".bin"

		 	# create binary file where each sample is a 32-bit float
			#fm_demod.astype('float32').tofile(fm_demod_fname)

		 	# plot PSD of selected block after extracting mono audio
		 	# ... change as needed

		 	# plot PSD of selected block after downsampling mono audio
		 	# ... change as needed

		 	# save figure to file
		 	#fig.savefig("../data/fmMonoBlock" + str(block_count) + ".png")

		block_count += 1

	print('Finished processing all the blocks from the recorded I/Q samples')

	# write audio data to file
	out_fname = r"data\fmMonoBlock.wav"
	wavfile.write(out_fname, int(audio_Fs), np.int16((audio_data/2)*32767))
	print("Written audio samples to \"" + out_fname + "\" in signed 16-bit format")

	# uncomment assuming you wish to show some plots
	# plt.show()
