/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Copyright by Nicola Nicolici
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "../include/filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"

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

	int rf_fs;
	int rf_fc = 100000;
	int rf_taps = 151;
	int rf_decim = 10;

	int audio_fs;
	int audio_decim = 5;
	int audio_taps = 101;
	int audio_fc = 16000;

	if (mode == 0){
		rf_fs = 2400000;
		audio_fs = 48000;
		
		// read from .raw

		std::vector<float> rf_coeff;
		impulseResponseLPF(rf_fs, rf_fc, rf_taps, rf_coeff);

		int block_size = 1024 * rf_decim * audio_decim * 2;
		int block_count = 0;

		std::vector<float> state_i_lpf_100k(rf_taps-1, 0.0);
		std::vector<float> state_q_lpf_100k(rf_taps-1, 0.0);
		std::vector<float> i_filt, q_filt;

		float prevI = 0;
		float prevQ = 0;
		std::vector<float> audio_data;
		std::vector<float> i_data_slice, q_data_slice;

		for (;; block_count++){
			std::vector<float> iq_block(block_size);
			readStdinBlockData(block_size, block_count, iq_block);
			if((std::cin.rdstate()) != 0){
				std::cerr << "End of input stream reached!" << std::endl;
				exit(1);
			}
			std::cerr << "Read block " << block_count << std::endl;

			auto start = iq_block.begin() + (block_count*block_size);
			auto end = iq_block.begin() + ((block_count+1)*block_size);
			copy(start, end, i_data_slice.begin());
			
			convolveFIR(i_filt, i_data_slice, state_i_lpf_100k);
			convolveFIR(q_filt, q_data_slice, state_q_lpf_100k);

			// decim

			demodFM(i_ds, q_ds, fm_demod, prevI, prevQ);

		}
	}
	else if (mode == 1){
		rf_fs = 1152000;
		audio_fs = 48000;

	}
	else if (mode == 2){
		rf_fs = 2400000;
		audio_fs = 44100;

	}
	else if (mode == 3){
		rf_fs = 2304000;
		audio_fs = 44100;
		
	}

	return 0;
}
