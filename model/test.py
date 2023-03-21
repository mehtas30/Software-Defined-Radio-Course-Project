import matplotlib.pyplot as plt
import numpy as np
import math

def lp_impulse_response_coeff(Fc, Fs, N_taps):
	normCutoff = Fc/(Fs/2)

	coefficients = np.zeros(N_taps)
	
	#print(f'{Fc}, {Fs}, {N_taps}')

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

def lp_filter(coefficients, data):

	coeff_len = len(coefficients)
	data_len = len(data)
	filtered_data = np.zeros(data_len)

	# discrete convolution
	for n in range(data_len):
		for k in range(coeff_len):
			if n-k >= 0 and n-k < data_len:
				filtered_data[n] += coefficients[k] * data[n-k]
	
	return filtered_data

def generateSin(Fs, interval, frequency = 7.0, amplitude = 5.0, phase = 0.0):
      dt = 1.0/Fs
      time = np.arange(0, interval, dt)
      x = amplitude * np.sin(2 * math.pi * frequency * time + phase)
      return time, x

f = 5.0
fs = 100.0
interval = 1.0
fc = 10.0
taps = 30

t, x = generateSin(fs, interval, frequency=f, amplitude=1.0, phase=0.0)
og = plt
og.plot(t, x)

x_pad = np.zeros(len(x) * 3)
for i in range(len(x)):
    x_pad[i*3] = x[i]
t_pad = np.arange(0, 3 * interval, 1/fs)	
padded = plt
padded.plot(t_pad, x_pad)

coeff = lp_impulse_response_coeff(fc, fs, taps)
x_us = 3*lp_filter(coeff, x_pad)
us = plt
us.plot(t_pad, x_us)

og.show()
padded.show()
us.show()