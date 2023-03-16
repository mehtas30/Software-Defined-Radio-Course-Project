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

void fmMonoProcessing(	int rf_fs, 	  int rf_fc,    int rf_taps,    int rf_decim, 
						int audio_fs, int audio_fc, int audio_taps, int audio_decim, 
						std::vector<float> &processed_data,
						int block_size, int block_count) 
{

	// coefficients for IQ -> IF LPFs, Fc = 100kHz
	std::vector<float> rf_coeff;
	impulseResponseLPF(rf_fs, rf_fc, rf_taps, rf_coeff);
	
	

	// coefficients for IF -> audio LPF, Fc = 16kHz
	std::vector<float> audio_coeff;
	impulseResponseLPF(240000, audio_fc, audio_taps, audio_coeff);
	

	// state saving for extracting FM band (Fc = 100kHz)
	std::vector<float> state_i_lpf_100k(rf_taps-1, 0.0);
	std::vector<float> state_q_lpf_100k(rf_taps-1, 0.0);
	// state saving for FM demodulation
	float prevI = 0;
	float prevQ = 0;
	// state saving for extracting mono audio (Fc = 16kHz)
	std::vector<float> audio_state(audio_taps-1, 0.0);

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
		for (unsigned int i = 0; i < iq_block.size(); i+=2){
			i_block[j] = iq_block[i];
			q_block[j] = iq_block[i+1];
			j++;
		};
	
		
		// LPF (Fc = 100kHz) extract FM band
		LPFilter(i_filt, i_block, rf_coeff, state_i_lpf_100k);
		LPFilter(q_filt, q_block, rf_coeff, state_q_lpf_100k);
		
		

		// from 2.4MS/s -> 240kS/s (decim=10)
		std::vector<float> i_ds;
		std::vector<float> q_ds;
		downSample(i_filt, i_ds, rf_decim);
		downSample(q_filt, q_ds, rf_decim);
		if (block_count == 0){
			for (int i = 0; i<30; i++){
				std::cerr << i_ds[i] << "," << q_ds[i] << "\n";
			}
		}
			
		// FM demodulation
		std::vector<float> fm_demod;
		demodFM(i_ds, q_ds, fm_demod, prevI, prevQ);

		// LPF (Fc = 16kHz)
		std::vector<float> audio_filt;
		LPFilter(audio_filt, fm_demod, audio_coeff, audio_state);

		// from 240kS/s -> 48kS/s
		processed_data.clear();
		downSample(audio_filt, processed_data, audio_decim);
		
		
		
		if (block_count == 10){
			std::vector<float> vector_index;
			genIndexVector(vector_index, fm_demod.size());
			//logVector("demod_time", vector_index, fm_demod);
			//logVector("i_block_time", vector_index, i_block);
			//logVector("q_block_time", vector_index, q_block);
			//logVector("iq_block_time", vector_index, iq_block);
			logVector("i_filt_time", vector_index, i_filt);
			logVector("q_filt_time", vector_index, q_filt);
			logVector("audio_filt_time", vector_index, audio_filt);
			
			std::vector<std::complex<float>> Xf;
			DFT(fm_demod, Xf);
			
			std::vector<float> Xmag;
			computeVectorMagnitude(Xf, Xmag);
			
			vector_index.clear();
			genIndexVector(vector_index, Xmag.size());
			logVector("demod_freq", vector_index, Xmag);
			
			std::vector<float> freq, psd_est;
			estimatePSD(freq, psd_est, fm_demod, block_size, 240);
			logVector("demod_psd", freq, psd_est);
		}
		
		// writing by block to stdout
		std::vector<short int> audio_data(block_size / rf_decim / audio_decim);
		for (unsigned int k=0; k < processed_data.size(); k++){
			if (std::isnan(processed_data[k])) audio_data[k] = 0;
			// prepare a block of audio data to be redirected to stdout at once
			else audio_data[k] = static_cast<short int>(processed_data[k] * 16384);
		}
		
		fwrite(&audio_data[0], sizeof(short int), audio_data.size(), stdout);
		//std::cerr << "Write block " << block_count << std::endl;
	}
}

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

	int audio_fs = 48000;
	int audio_decim = 5;
	int audio_taps = 101;
	int audio_fc = 16000;
	std::vector<float> processed_data;

	int block_size = 1024 * rf_decim * audio_decim * 2;
	int block_count = 0;

	if (mode == 0){ 	 rf_fs = 2400000; audio_fs = 48000; }
	else if (mode == 1){ rf_fs = 1152000; audio_fs = 48000; }
	else if (mode == 2){ rf_fs = 2400000; audio_fs = 44100; }
	else if (mode == 3){ rf_fs = 2304000; audio_fs = 44100; }

	fmMonoProcessing(rf_fs, rf_fc, rf_taps, rf_decim, audio_fs, audio_fc, audio_taps, audio_decim, processed_data, block_size, block_count);

	std::cerr << "done" << std::endl;
	return 0;
}

// cat ../data/samples_u8.raw | ./project | aplay -c 1 -f S16_LE -r 48000
// gnuplot -e 'set terminal png size 1024,768' ../data/example.gnuplot > ../data/example.png
