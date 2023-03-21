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
	
	// extracted 0-100kHz band
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
		for (int i = 0; i < (int)iq_block.size(); i+=2){
			i_block[j] = iq_block[i];
			q_block[j] = iq_block[i+1];
			j++;
		};
		
		// LPF (Fc = 100kHz) extract FM band
		LPFilter(i_filt, state_i_lpf_100k, i_block, rf_coeff);
		LPFilter(q_filt, state_q_lpf_100k, q_block, rf_coeff);	
		
		
		// from 2.4MS/s -> 240kS/s (decim=10)
		std::vector<float> i_ds;
		std::vector<float> q_ds;
		downsample(i_ds, i_filt, rf_decim);
		downsample(q_ds, q_filt, rf_decim);
		
		//if (block_count == 0){
			//std::cerr << "\n";
			//for (int i = 0; i < 30; i++){
				//std::cerr << "i=" << i_block[i] << ", q=" << q_block[i] << "\n";
			//}
			//std::cerr << "\n" << std::endl;
			//for (int i = 0; i < 30; i++){
				//std::cerr << "i_filt=" << i_filt[i] << ", q_filt=" << q_filt[i] << "\n";
			//}
			//std::cerr << "\n" << std::endl;
			//for (int i = 0; i < 30; i++){
				//std::cerr << "i_ds=" << i_ds[i] << ", q_ds=" << q_ds[i] << "\n";
			//}
			//std::cerr << "\n" << std::endl;
		//}
			
		// FM demodulation
		std::vector<float> fm_demod;
		FMDemod(fm_demod, prev_i, prev_q, i_ds, q_ds);
		
		// zero padding for upsampling
		std::vector<float> fm_demod_us;
		upsample(fm_demod_us, fm_demod, audio_interp);

		// LPF (Fc = 16kHz)
		std::vector<float> audio_filt;
		LPFilter(audio_filt, audio_state, fm_demod_us, audio_coeff);

		// from 240kS/s -> 48kS/s
		downsample(processed_data, audio_filt, audio_decim);

		
		//if (block_count == 0){
			//std::cerr << "\n";
			//for (int i = 0; i < 30; i++){
				//std::cerr << "fm_demod=" << fm_demod[i] << "\n";
			//}
			//std::cerr << "\n" << std::endl;
			//for (int i = 0; i < 30; i++){
				//std::cerr << "audio_filt=" << audio_filt[i] << "\n";
			//}
			//std::cerr << "\n" << std::endl;
			//for (int i = 0; i < 30; i++){
				//std::cerr << "processed_data=" << processed_data[i] << "\n";
			//}
			//std::cerr << "\n" << std::endl;
		//}
		
		//if (block_count == 10){
			//std::vector<float> vector_index;
			//genIndexVector(vector_index, fm_demod.size());
			//logVector("demod_time", vector_index, fm_demod);
			//logVector("i_block_time", vector_index, i_block);
			//logVector("q_block_time", vector_index, q_block);
			//logVector("iq_block_time", vector_index, iq_block);
			//logVector("i_filt_time", vector_index, i_filt);
			//logVector("q_filt_time", vector_index, q_filt);
			//logVector("audio_filt_time", vector_index, audio_filt);
			
			//std::vector<std::complex<float>> Xf;
			//DFT(fm_demod, Xf);
			
			//std::vector<float> Xmag;
			//computeVectorMagnitude(Xf, Xmag);
			
			//vector_index.clear();
			//genIndexVector(vector_index, Xmag.size());
			//logVector("demod_freq", vector_index, Xmag);
			
			//std::vector<float> freq, psd_est;
			//estimatePSD(freq, psd_est, fm_demod, block_size, 240);
			//logVector("demod_psd", freq, psd_est);
		//}
		
		
		// writing by block to stdout
		std::vector<short int> audio_data(block_size / rf_decim / audio_decim);
		for (unsigned int k=0; k < processed_data.size(); k++){
			if (std::isnan(processed_data[k])) audio_data[k] = 0;
			// prepare a block of audio data to be redirected to stdout at once
			else audio_data[k] = static_cast<short int>(processed_data[k] * 16384);
		}
		
		fwrite(&audio_data[0], sizeof(short int), audio_data.size(), stdout);
	}
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
	int rf_taps = 30;
	int rf_decim = 10;
	
	int if_fs = 240000;

	int audio_fs = 48000;
	int audio_fc = 16000;
	int audio_taps = 30;
	int audio_decim = 5;
	int audio_interp = 1;
	
	std::vector<float> processed_data;

	int block_size = 128 * rf_decim * 5 * 2;
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
