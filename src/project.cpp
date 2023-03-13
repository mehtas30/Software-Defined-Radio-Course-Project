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

int main(int argc, char* argv[])
{
	// binary files can be generated through the
	// Python models from the "../model/" sub-folder
	const std::string in_fname = "../data/fm_demod_10.bin";
	std::vector<float> bin_data;
	readBinData(in_fname, bin_data);

	int mode = 0;

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
	int rf_taps = 151;
	int rf_decim = 10;

	int audio_fs;
	int audio_decim = 5;
	int audio_taps = 101;
	int audio_fc = 16000;
	std::vector<float> audio_data;

	int block_size = 1024 * rf_decim * audio_decim * 2;
	int block_count = 0;

	if (mode == 0){ 	 rf_fs = 2400000; audio_fs = 48000; }
	else if (mode == 1){ rf_fs = 1152000; audio_fs = 48000; }
	else if (mode == 2){ rf_fs = 2400000; audio_fs = 44100; }
	else if (mode == 3){ rf_fs = 2304000; audio_fs = 44100; }

	// read from .raw

	fmMonoProcessing(rf_fs, rf_fc, rf_taps, rf_decim, audio_fs, audio_decim, audio_taps, audio_fc, audio_data);

	std::vector<short int> final_audio(block_size);
	for (int i=0; i<audio_data.size(); i++) {
		if (std::isnan(audio_data[i])) final_audio[i] = 0;
		else final_audio[i] = static_cast<short int>(audio_data[i] * 16384);
	}

	fwrite(&final_audio[0], sizeof(short int), final_audio.size(), stdout);

	return 0;
}

void fmMonoProcessing(int rf_fs, int rf_fc, int rf_taps, int rf_decim, int audio_fs, int audio_decim, int audio_taps, int audio_fc, std::vector<float> &audio_data) {
	int block_size = 1024 * rf_decim * audio_decim * 2;
	int block_count = 0;

	// coefficients for IQ -> IF LPFs, Fc = 100kHz
	std::vector<float> rf_coeff;
	impulseResponseLPF(rf_fs, rf_fc, rf_taps, rf_coeff);

	// coefficients for IF -> audio LPF, Fc = 16kHz
	std::vector<float> audio_coeff;
	impulseResponseLPF(audio_fs, audio_fc, audio_taps, audio_coeff);

	// state saving for extracting FM band (Fc = 100kHz)
	std::vector<float> state_i_lpf_100k(rf_taps-1, 0.0);
	std::vector<float> state_q_lpf_100k(rf_taps-1, 0.0);
	// state saving for FM demodulation
	float prevI = 0;
	float prevQ = 0;
	// state saving for extracting mono audio (Fc = 16kHz)
	std::vector<float> audio_state(audio_taps-1, 0.0);

	// audio buffer that stores all audio blocks

	std::vector<float> i_filt;
	std::vector<float> q_filt;

	for (;; block_count++){
		std::vector<float> iq_block(block_size * 2);
		readStdinBlockData(block_size * 2, block_count, iq_block);
		if((std::cin.rdstate()) != 0){
			std::cerr << "End of input stream reached!" << std::endl;
			exit(1);
		}
		std::cerr << "Read block " << block_count << std::endl;

		// separate iq_block into I and Q
		std::vector<float> i_block(block_size);
		std::vector<float> q_block(block_size);
		int j = 0;
		for (int i = 0; i < iq_block.size(); i+=2){
			i_block[j] = iq_block[i];
			q_block[j] = iq_block[i+1];
			j++;
		}
		
		LPFilter(i_filt, i_block, rf_coeff, state_i_lpf_100k);
		LPFilter(q_filt, q_block, rf_coeff, state_q_lpf_100k);

		// decim
		std::vector<float> i_ds;
		std::vector<float> q_ds;

		i_ds = downSample(i_filt, rf_decim);
		q_ds = downSample(q_filt, rf_decim);

		std::vector<float> fm_demod;
		demodFM(i_ds, q_ds, fm_demod, prevI, prevQ);

		std::vector<float> audio_filt;
		LPFilter(audio_filt, fm_demod, audio_coeff, audio_state);

		std::vector<float> audio_block;
		audio_block = downSample(audio_filt,audio_decim);

		audio_data.insert(audio_data.end(), audio_block.begin(), audio_block.end());
	}
}

std::vector<float> downSample(std::vector<float> original, int decim) {
	std::vector<float> downSample;

	for (int i = 0; i < original.size(); i += decim) {
		downSample.push_back(original[i]);
	}

	return downSample;
}