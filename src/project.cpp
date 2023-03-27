/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Copyright by Nicola Nicolici
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "../include/dy4.h"
#include "../include/filter.h"
#include "../include/fourier.h"
#include "../include/genfunc.h"
#include "../include/iofunc.h"
#include "../include/logfunc.h"
#include <chrono>

void fmMonoProcessing(	int rf_fs, 	  int rf_fc,    int rf_taps,    int rf_decim,	 int if_fs,
						int audio_fs, int audio_fc, int audio_taps, int audio_decim, int audio_interp,
						std::vector<float> &processed_data,
						int block_size, int block_count) 
{
	std::cerr << "rf_fs=" << rf_fs << "\nrf_fc=" << rf_fc << "\nrf_taps=" << rf_taps << "\nrf_decim=" << rf_decim << "\nif_fs=" << if_fs << "\naudio_fs=" << audio_fs << "\naudio_fc=" << audio_fc << "\naudio_taps=" << audio_taps << "\naudio_decim=" << audio_decim << "\naudio_interp=" << audio_interp << "\nblock_size=" << block_size << "\nblock_count=" << block_count << std::endl;

	// coefficients for IQ -> IF LPFs, Fc = 100kHz
	std::vector<float> rf_coeff;
	impulseResponseLPF(rf_coeff, rf_fs, rf_fc, rf_taps, 1);
	
	// coefficients for IF -> audio LPF, Fc = 16kHz
	std::vector<float> audio_coeff;
	impulseResponseLPF(audio_coeff, if_fs, audio_fc, audio_taps, audio_interp);
	
	// state saving for extracting FM band (Fc = 100kHz)
	std::vector<float> state_i_lpf_100k(rf_taps-1, 0.0);
	std::vector<float> state_q_lpf_100k(rf_taps-1, 0.0);
	// state saving for FM demodulation
	float prev_i = 0.0;
	float prev_q = 0.0;
	// state saving for extracting mono audio (Fc = 16kHz)
	std::vector<float> audio_state(audio_taps-1, 0.0);
	
	// split IQ block into I and Q
	std::vector<float> i_block;
	std::vector<float> q_block;
	i_block.reserve(block_size); i_block.resize(block_size);
	q_block.reserve(block_size); q_block.resize(block_size);

	// downsampled I and Q
	std::vector<float> i_ds; 
	std::vector<float> q_ds; 
	
	// demodulated data
	std::vector<float> fm_demod; 
	
	for (;; block_count++){
		
		std::vector<float> iq_block(block_size * 2);
		readStdinBlockData(block_size * 2, block_count, iq_block);
		if((std::cin.rdstate()) != 0){
			std::cerr << "End of input stream reached!" << std::endl;
			exit(1);
		}
		std::cerr << "Read block " << block_count << std::endl;

		// separate iq_block into I and Q
		int j = 0;
		for (int i = 0; i < (int)iq_block.size(); i+=2){
			i_block[j] = iq_block[i];
			q_block[j] = iq_block[i+1];
			j++;
		};

		// LPF (Fc = 100kHz) extract FM band
		//std::cerr << "Extracting FM band...\ntaps=" << rf_coeff.size() << ", block size=" << i_block.size() << std::endl;
		//LPFilter(i_filt, state_i_lpf_100k, i_block, rf_coeff);
		//LPFilter(q_filt, state_q_lpf_100k, q_block, rf_coeff);	
		
		// from 2.4MS/s -> 240kS/s (decim=10)
		//std::cerr << "Downsampling FM band..." << std::endl;
		//downsample(i_ds, i_filt, rf_decim);
		//downsample(q_ds, q_filt, rf_decim);
		
		// resampling filter (2.4MS/s -> 240kS/s, 0 - 100kHz)
		resample(i_ds, state_i_lpf_100k, i_block, rf_coeff, 1, rf_decim);
		resample(q_ds, state_q_lpf_100k, q_block, rf_coeff, 1, rf_decim);
		
		// FM demodulation
		FMDemod(fm_demod, prev_i, prev_q, i_ds, q_ds);
		
		// resampling filter (240kS/s -> 48kS/s, 0 - 16kHz)
		resample(processed_data, audio_state, fm_demod, audio_coeff, audio_interp, audio_decim);
		
		// zero padding for upsampling
		//std::cerr << "Upsampling..." << std::endl;
		//upsample(fm_demod_us, fm_demod, audio_interp);

		// LPF (Fc = 16kHz)
		//std::cerr << "Extracting mono audio...\ntaps=" << audio_coeff.size() << ", block size=" << fm_demod_us.size() << std::endl;
		//LPFilter(audio_filt, audio_state, fm_demod_us, audio_coeff);
		
		// from 240kS/s -> 48kS/s
		//std::cerr << "Downsampling mono audio..." << std::endl;
		//downsample(processed_data, audio_filt, audio_decim);
		
		// writing by block to stdout
		//std::cerr << "Writing to stdout..." << std::endl;
		std::vector<short int> audio_data;
		audio_data.reserve(processed_data.size());
		audio_data.resize(processed_data.size(), 0.0);
		
		for (unsigned int k=0; k < processed_data.size(); k++){
			if (std::isnan(processed_data[k])) audio_data[k] = 0;
			// prepare a block of audio data to be redirected to stdout at once
			else audio_data[k] = static_cast<short int>(processed_data[k] * 16384);
		}
		
		fwrite(&audio_data[0], sizeof(short int), audio_data.size(), stdout);
	}
}
void stereo(int rf_fs, int rf_fc, int rf_taps, int rf_decim, int if_fs,
			int audio_fs, int audio_fc, int audio_taps, int audio_decim, int audio_interp,
			std::vector<float> &processed_data,
			int block_size, int block_count)
{
	std::cerr << "rf_fs=" << rf_fs << "\nrf_fc=" << rf_fc << "\nrf_taps=" << rf_taps << "\nrf_decim=" << rf_decim << "\nif_fs=" << if_fs << "\naudio_fs=" << audio_fs << "\naudio_fc=" << audio_fc << "\naudio_taps=" << audio_taps << "\naudio_decim=" << audio_decim << "\naudio_interp=" << audio_interp << "\nblock_size=" << block_size << "\nblock_count=" << block_count << std::endl;
	// coefficients for IQ -> IF LPFs, Fc = 100kHz
	std::vector<float> rf_coeff;
	impulseResponseLPF(rf_coeff, rf_fs, rf_fc, rf_taps, 1);
	// coefficients for IF -> audio LPF, Fc = 16kHz
	std::vector<float> mono_coeff;
	impulseResponseLPF(mono_coeff, 240000, audio_fc, audio_taps, 1);
	// state saving for extracting FM band, Fc=100Khz
	std::vector<float> state_i_lpf_100k(rf_taps - 1, 0.0);
	std::vector<float> state_q_lpf_100k(rf_taps - 1, 0.0);
	// state saving for FM demod
	float prev_i = 0.0;
	float prev_q = 0.0;
	// extracting mono state saving
	std::vector<float> monoState(audio_taps - 1, 0.0);
	std::vector<float> channelExtractState(audio_taps - 1, 0.0);
	std::vector<float> carrierRecoveryState(audio_taps - 1, 0.0);
	// audio buffer that stores all the audio blocks
	std::vector<float> carrierRecoveryData;
	std::vector<float> channelExtractData;
	std::vector<float> monoData;
	std::vector<float> leftData;
	std::vector<float> rightData;
	// split IQ block into I and Q
	std::vector<float> i_block;
	std::vector<float> q_block;
	i_block.reserve(block_size);
	i_block.resize(block_size);
	q_block.reserve(block_size);
	q_block.resize(block_size);

	// downsampled I and Q
	std::vector<float> i_ds;
	std::vector<float> q_ds;

	// demodulated data
	std::vector<float> fm_demod;
	// extract coeffs and state
	std::vector<float> channelExtractCoeff;
	// recovery
	std::vector<float> carrierRecovCoeff;
	// process
	std::vector<float> mixerdata;
	std::vector<float> downsampledStereo;
	std::vector<float> upsampledStereo;
	std::vector<float> leftNewData;
	std::vector<float> rightNewData;
	for (;; block_count++)
	{
		// input
		std::vector<float> iq_block(block_size * 2);
		readStdinBlockData(block_size * 2, block_count, iq_block);
		if ((std::cin.rdstate()) != 0)
		{
			std::cerr << "End of input stream reached!" << std::endl;
			exit(1);
		}
		std::cerr << "Read block " << block_count << std::endl;

		// separate iq_block into I and Q
		int j = 0;
		for (int i = 0; i < (int)iq_block.size(); i += 2)
		{
			i_block[j] = iq_block[i];
			q_block[j] = iq_block[i + 1];
			j++;
		};
		// resampling filter (2.4MS/s -> 240kS/s, 0 - 100kHz)
		resample(i_ds, state_i_lpf_100k, i_block, rf_coeff, 1, rf_decim);
		resample(q_ds, state_q_lpf_100k, q_block, rf_coeff, 1, rf_decim);
		// Fm demod
		FMDemod(fm_demod, prev_i, prev_q, i_ds, q_ds);
		// stereo extract
		bandpassfilter(channelExtractCoeff, 22000.0, 45000.0, 240000.0, audio_taps); // why error
		LPFilter(channelExtractData, channelExtractState, fm_demod, channelExtractCoeff);
		resample(channelExtractData, channelExtractState, fm_demod, channelExtractCoeff, 1, 1);
		// recovery
		bandpassfilter(carrierRecovCoeff, 18500, 19500, 240000, audio_taps);
		LPFilter(carrierRecoveryData, carrierRecoveryState, fm_demod, carrierRecovCoeff);
		fmPLL(carrierRecoveryData, 19000, 240000, 2, 0, 0.01);
		resample(carrierRecoveryData, carrierRecoveryState, fm_demod, carrierRecovCoeff, 1, 1);
		// mono
		resample(monoData, monoState, fm_demod, mono_coeff, audio_interp, audio_decim);
		// process
		mixer(mixerdata, channelExtractData, carrierRecoveryData);
		downsample(downsampledStereo, mixerdata, audio_decim);
		upsample(upsampledStereo, downsampledStereo, audio_interp);
		lrExtraction(leftNewData, rightNewData, monoData, upsampledStereo);
		leftData.insert(leftData.end(), leftNewData.begin(), leftNewData.end());
		rightData.insert(rightData.end(), rightNewData.begin(), rightNewData.end());
	}
	// writing by block to stdout
	// std::cerr << "Writing to stdout..." << std::endl;
	std::vector<short int> audioleft_data;
	audioleft_data.reserve(leftData.size());
	audioleft_data.resize(leftData.size(), 0.0);
	std::vector<short int> audioright_data;
	audioright_data.reserve(rightData.size());
	audioright_data.resize(rightData.size(), 0.0);
	for (unsigned int k = 0; k < leftData.size(); k++)
	{
		if (std::isnan(leftData[k]))
			audioleft_data[k] = 0;
		// prepare a block of audio data to be redirected to stdout at once
		else
			audioleft_data[k] = static_cast<short int>(leftData[k] * 16384);
	}

	fwrite(&audioright_data[0], sizeof(short int), audioright_data.size(), stdout);

	for (unsigned int k = 0; k < rightData.size(); k++)
	{
		if (std::isnan(rightData[k]))
			audioright_data[k] = 0;
		// prepare a block of audio data to be redirected to stdout at once
		else
			audioright_data[k] = static_cast<short int>(rightData[k] * 16384);
	}

	fwrite(&audioright_data[0], sizeof(short int), audioright_data.size(), stdout);
}
int main(int argc, char* argv[])
{
	// binary files can be generated through the
	// Python models from the "../model/" sub-folder
	const std::string in_fname = "../data/fm_demod_10.bin";
	std::vector<float> bin_data;
	readBinData(in_fname, bin_data);

	int mode = 2;

	if (argc < 2){
		std::cerr << "Operating in default mode 0" << std::endl;
	}
	else if (argc == 2){
		mode = atoi(argv[1]);
		if (mode > 3){
			std::cerr << "Wrong mode " << mode << std::endl;
			exit(1);
		}
	}
	else{
		std::cerr << "Usage: " << argv[0] << std::endl;
		std::cerr << "or " << std::endl;
		std::cerr << "Usage " << argv[0] << " <mode>" << std::endl;
		std::cerr << "\t\t <mode> is a value from 0 to 3" << std::endl;
		exit(1);
	}

	std::cerr << "Operating in mode " << mode << std::endl;

	int rf_fs = 2400000;
	int rf_fc = 100000;
	int rf_taps = 51;
	int rf_decim = 10;
	
	int if_fs = 240000;

	int audio_fs = 48000;
	int audio_fc = 16000;
	int audio_taps = 31;
	int audio_decim = 5;
	int audio_interp = 1;
	
	std::vector<float> processed_data;

	int block_size = 128 * rf_decim * audio_decim;
	int block_count = 0;

	if (mode == 0){
		rf_fs = 2400000; rf_decim = 10;
		if_fs = 240000;
		audio_fs = 48000; audio_decim = 5;
	}
	else if (mode == 1){ 
		rf_fs = 1152000; rf_decim = 4;
		if_fs = 288000;
		audio_fs = 48000; audio_decim = 4;
	}
	else if (mode == 2){ 
		rf_fs = 2400000; rf_decim = 10;
		if_fs = 240000;
		audio_fs = 44100; audio_decim = 800; audio_interp = 147;
		audio_taps *= audio_interp;
		if_fs *= audio_interp;
	}
	else if (mode == 3){ 
		rf_fs = 2304000; rf_decim = 9;
		if_fs = 256000;
		audio_fs = 44100; audio_decim = 2560; audio_interp = 441;
		audio_taps *= audio_interp;
		if_fs *= audio_interp;
	}

	fmMonoProcessing(rf_fs, rf_fc, rf_taps, rf_decim, if_fs, 
					audio_fs, audio_fc, audio_taps, audio_decim, audio_interp, 
					processed_data, block_size, block_count);

	std::cerr << "done" << std::endl;
	return 0;
}

// cat ../data/samples_u8.raw | ./project | aplay -c 1 -f S16_LE -r 48000
// gnuplot -e 'set terminal png size 1024,768' ../data/example.gnuplot > ../data/example.png
