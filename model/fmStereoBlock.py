import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy import signal
import numpy as np
import math

from fmSupportLib import fmDemodArctan, fmPlotPSD

def bandpassFilt(fb, fe, fs, n_taps):
    normCent = ((fe + fb)/2) / (fs/2)
    normPass = (fe - fb) / (fs/2)

    coefficients = np.zeros(n_taps-1)

    for i in range (n_taps-1):
        if (i == (n_taps-1)/2):
            coefficients[i] = normPass
        else:
            value = math.pi * (normPass/2) * (i - ((n_taps-1)/2))
            coefficients[i] = normPass * (math.sin(value) / value)

        coefficients[i] *= math.cos(i * math.pi * normCent)
        coefficients[i] *= math.sin((i * math.pi) / n_taps)**2

    return coefficients

def fmPll(pllIn, freq, Fs, nocoScale = 1.0, phaseAdjust = 0.0, normBandwidth = 0.01):


    Cp = 2.666
    Ci = 3.555

    Kp = (normBandwidth) * Cp
    Ki = (normBandwidth*normBandwidth) * Ci

    ncoOut = np.empty(len(pllIn) + 1)

    integrator = 0.0
    phaseEst = 0.0
    feedbackI = 1.0
    feedbackQ = 0.0
    ncoOut[0] = 1.0
    trigOffset = 0

    for i in range (len(pllIn)):
        errorI = pllIn[i] * (+feedbackI)
        errorQ = pllIn[i] * (-feedbackQ)

        errorD = math.atan2(errorQ, errorI)

        integrator = integrator + Ki*errorD

        phaseEst = phaseEst + Kp*errorD + integrator

        trigOffset += 1
        trigArg = 2*math.pi * (freq/Fs) * (trigOffset) + phaseEst
        feedbackI = math.cos(trigArg)
        feedbackQ = math.sin(trigArg)
        ncoOut[i+1] = math.cos(trigArg * nocoScale + phaseAdjust)

    return ncoOut

def filter(coefficients, data, state):

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

def lp_impulse_response_coeff(Fc, Fs, N_taps):
    normCutoff = Fc/(Fs/2)

    coefficients = np.zeros(N_taps)

    for i in range(N_taps):
        if (i == ((N_taps-1)/2)):
            coefficients[i] = normCutoff

        else:
            # sinc = sin(x)/x
            numerator = math.sin(math.pi * normCutoff * (i - ((N_taps-1) / 2)))
            denominator = math.pi * normCutoff * (i - ((N_taps-1) / 2))

            # normalize coefficients
            coefficients[i] = normCutoff * (numerator/denominator)

        # window the coefficients
        coefficients[i] = coefficients[i] * ((math.sin((i*math.pi) / N_taps)) ** 2)

    return coefficients

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

def mixer(arr1, arr2):
    mixedData = []
    for i in range(len(arr2)):
            mixedData[i] = arr1[i] * arr2[i]
    return mixedData

def lrExtraction(monoData, stereoData):
    leftData = []
    rightData = []
    for i in range(len(monoData)):
        leftData = (monoData[i] + stereoData[i]) / 2
        rightData = (monoData[i] - stereoData[i]) / -2
        
rf_Fs = 2.4e6
rf_Fc = 100e3
rf_taps = 151
rf_decim = 10

audio_Fs = 48e3
audio_Fc = 16e3
audio_taps = 101
audio_decim = 5

if __name__ == "__main__":

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

    # select a block_size that is a multiple of KB
    # and a multiple of decimation factors
    block_size = 1024 * rf_decim * audio_decim * 2
    block_count = 0

    # coefficients for IQ -> IF LPFs, Fc = 100kHz
    rf_coeff = lp_impulse_response_coeff(rf_Fc, rf_Fs, rf_taps)

    # coefficients for IF -> audio LPF, Fc = 16kHz
    audio_coeff = lp_impulse_response_coeff(audio_Fc, 240e3, audio_taps)

    # state-saving for extracting FM band (Fc = 100kHz)
    state_i_lpf_100k = np.zeros(rf_taps-1)
    state_q_lpf_100k = np.zeros(rf_taps-1)
    # state-saving for FM demodulation
    prevI = 0
    prevQ = 0
    # state-saving for extracting mono audio (Fc = 16kHz)
    channelExtractState = np.zeros(audio_taps-1)
    carrierRecoveryState = np.zeros(audio_taps-1)

    # audio buffer that stores all the audio blocks
    carrierRecoveryData = np.array([]) # used to concatenate filtered blocks (audio data)
    channelExtractData = np.array([])

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
        #i_filt, state_i_lpf_100k = lp_filter(rf_coeff,
        #		       iq_data[(block_count)*block_size+1:(block_count+1)*block_size:2],
        #			   state_i_lpf_100k)
        #q_filt, state_q_lpf_100k = lp_filter(rf_coeff,
        #		       iq_data[(block_count)*block_size+1:(block_count+1)*block_size:2],
        #			   state_q_lpf_100k)


        # downsample the I/Q data from the FM channel
        i_ds = i_filt[::rf_decim]
        q_ds = q_filt[::rf_decim]

        # FM demodulator
        # you will need to implement your own FM demodulation based on:
        # https://www.embedded.com/dsp-tricks-frequency-demodulation-algorithms/
        # see more comments on fmSupportLib.py - take particular notice that
        # you MUST have also "custom" state-saving for your own FM demodulator
        fm_demod, prevI, prevQ = myDemod(i_ds, q_ds, prevI, prevQ)

        """
        Stereo Carrier Recovery
        """
        carrierRecoveryCoeff = bandpassFilt(18.5e3, 19.5e3, 240e3, audio_taps)

        carrierRecoveryFiltered, carrierRecoveryState = filter(carrierRecoveryCoeff, fm_demod, carrierRecoveryState)

        carrierRecoveryBlock = fmPll(carrierRecoveryFiltered,19e3,240e3,2,0,0.01)

        if (block_count != 0):
            carrierRecoveryData = np.concatenate([carrierRecoveryData, carrierRecoveryBlock])
        else:
            carrierRecoveryData = carrierRecoveryBlock

        """
        Stereo Channel Extraction
        """

        channelExtractCoeff = bandpassFilt(22e3, 54e3, 240e3, audio_taps)

        channelExtractFiltered, channelExtractState = filter(channelExtractCoeff, fm_demod, channelExtractState)

        if (block_count == 0):
            channelExtractData = channelExtractFiltered
        else:
<<<<<<< HEAD
            channelExtractData = np.concatenate(channelExtractData, channelExtractFiltered)
        """
        
        #Stereo Processing
        #Mixer
        mixedData = mixer(channelExtractData, carrierRecoveryData)
        #LPF
        low_pass_coeff = lp_impulse_response_coeff(audio_Fc, audio_Fs, audio_taps)
        stereoData = lp_filter(low_pass_coeff, mixedData)
        #LR blocks
        lData,rData = lrExtraction(monoData,stereoData)
        
=======
            channelExtractData = np.concatenate([channelExtractData, channelExtractFiltered])
>>>>>>> db0d12ba79385d4bef3d890fae7bb159860b80ae


        block_count += 1
