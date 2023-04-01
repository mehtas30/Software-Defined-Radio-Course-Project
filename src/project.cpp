/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Copyright by Nicola Nicolici
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "../include/dy4.h"
#include "../include/filter.h"
#include "../include/iofunc.h"

void monoProcessing(std::vector<float> &mono_data,
					const std::vector<float> &audio_coeff, std::vector<float> &audio_state,
					const std::vector<float> &demod_data,
					const int audio_decim, const int audio_interp, const int block_count) 
{
	// resampling filter (240kS/s -> 48kS/s, 0 - 16kHz)
	resample(mono_data, audio_state, demod_data, audio_coeff, audio_interp, audio_decim);

}


void stereoProcessing(std::vector<float> &stereo_data, 
					const std::vector<float> &carrier_coeff, std::vector<float> &carrier_state,
					const std::vector<float> &channel_coeff, std::vector<float> &channel_state,
					const std::vector<float> &audio_coeff, std::vector<float> &audio_state,
					const std::vector<float> &demod_data,
					const int audio_decim, const int audio_interp, const int if_fs, 
					float &integrator, float &phaseEst, float &feedbackI, float &feedbackQ, float &ncoOut_state, float &trigOffset,
					const int block_count)
{
	// channel extraction
	std::vector<float> channel_data;
	
	// carrier recovery
	std::vector<float> carrier_data;
	
	// stereo processing
	std::vector<float> mixer_data;

	// stereo channel extraction
	resample(channel_data, channel_state, demod_data, channel_coeff, 1, 1);

	// stereo carrier recovery
	resample(carrier_data, carrier_state, demod_data, carrier_coeff, 1, 1);
	PLL(carrier_data, 19000, if_fs, 2, 0, 0.01, integrator, phaseEst, feedbackI, feedbackQ, ncoOut_state, trigOffset);
	
	// stereo processing
	mixer(mixer_data, channel_data, carrier_data);
	
	// 16kHz filter
	resample(stereo_data, audio_state, mixer_data, audio_coeff, audio_interp, audio_decim);
}


int main(int argc, char* argv[])
{
	// binary files can be generated through the
	// Python models from the "../model/" sub-folder
	const std::string in_fname = "../data/fm_demod_10.bin";
	std::vector<float> bin_data;
	readBinData(in_fname, bin_data);

	int mode = 0;
	int channels = 1;

	if (argc < 3){
		std::cerr << "Operating in default mode 0, mono" << std::endl;
	}
	else if (argc == 3){
		mode = atoi(argv[1]);
		channels = atoi(argv[2]);
		if (mode < 0 || mode > 3){
			std::cerr << "Invalid mode: " << mode << "!" << std::endl;
			exit(1);
		}
		if (channels < 1 || channels > 2){
			std::cerr << "Invaild channel: " << channels << "!" << std::endl;
			exit(1);
		}
	}
	else{
		std::cerr << "Usage: " << argv[0] << std::endl;
		std::cerr << "or " << std::endl;
		std::cerr << "Usage " << argv[0] << " <mode> <channels>" << std::endl;
		std::cerr << "\t\t<mode> is a value from 0 to 3\n\t\t   <channels> is either 1 or 2" << std::endl;
		exit(1);
	}

	if (channels == 1) std::cerr << "Operating in mode " << mode << ", mono" << std::endl;
	else std::cerr << "Operating in mode " << mode << ", stereo" << std::endl;

	int rf_fs = 2400000;
	int rf_fc = 100000;
	int rf_taps = 51;
	int rf_decim = 10;
	
	int if_fs = 240000;
	int bp_fs = 240000;
	int bp_taps = 51;

	int audio_fc = 16000;
	int audio_taps = 51;
	int audio_decim = 5;
	int audio_interp = 1;
	
	int mono_delay = 5;

	if (mode == 0){
		rf_fs = 2400000; rf_decim = 10;
		if_fs = 240000;
		bp_fs = 240000;
		audio_decim = 5;
	}
	else if (mode == 1){ 
		rf_fs = 1152000; rf_decim = 4;
		if_fs = 288000;
		bp_fs = 288000;
		audio_decim = 6;
	}
	else if (mode == 2){ 
		rf_fs = 2400000; rf_decim = 10;
		if_fs = 240000;
		bp_fs = 240000;
		audio_decim = 800; audio_interp = 147;
		audio_taps *= audio_interp;
		if_fs *= audio_interp;
	}
	else if (mode == 3){ 
		rf_fs = 2304000; rf_decim = 9;
		if_fs = 256000;
		bp_fs = 256000;
		audio_decim = 2560; audio_interp = 441;
		audio_taps *= audio_interp;
		if_fs *= audio_interp;
	}
	
	int block_size = 256 * rf_decim * audio_decim;
	int block_count = 0;
	
	// begin processing
	
	// FM band extraction
	std::vector<float> iq_block(block_size);
	
	int i_block_size = int(iq_block.size() * 0.5);
	std::vector<float> i_block;
	std::vector<float> q_block;
	i_block.reserve(i_block_size); i_block.resize(i_block_size);
	q_block.reserve(i_block_size); q_block.resize(i_block_size);
	std::vector<float> state_i_lpf_100k(rf_taps-1, 0.0);
	std::vector<float> state_q_lpf_100k(rf_taps-1, 0.0);
	// coefficients for IQ -> IF LPFs, Fc = 100kHz
	std::vector<float> rf_coeff;
	impulseResponseLPF(rf_coeff, rf_fs, rf_fc, rf_taps, 1);
	
	// downsampled I and Q
	std::vector<float> i_ds; 
	std::vector<float> q_ds; 
	
	// FM demodulation
	float prev_i = 0.0;
	float prev_q = 0.0;
	std::vector<float> demod_data;
	
	// stereo channel extraction
	std::vector<float> channel_state(bp_taps-1, 0.0);
	// coefficients for stereo channel extraction, 22kHz to 54kHz
	std::vector<float> channel_coeff;
	impulseResponseBPF(channel_coeff, bp_fs, 22000.0, 54000.0, bp_taps);
	
	// stereo carrier recovery 
	std::vector<float> carrier_state(bp_taps-1, 0.0);
	// coefficients for stereo carrier recovery, 18.5kHz to 19.5kHz
	std::vector<float> carrier_coeff;
	impulseResponseBPF(carrier_coeff, bp_fs, 18500, 19500, bp_taps);
	// PLL state saving
	float integrator = 0.0;
	float phaseEst = 0.0;
	float feedbackI = 1.0;
	float feedbackQ = 0.0;
	float trigOffset = 0.0;
	float ncoOut_state = 1.0;
	
	// state saving for extracting audio (Fc = 16kHz)
	std::vector<float> audio_state(audio_taps-1, 0.0);
	// coefficients for IF -> audio LPF, Fc = 16kHz
	std::vector<float> audio_coeff;
	impulseResponseLPF(audio_coeff, if_fs, audio_fc, audio_taps, audio_interp);
	
	std::vector<float> mono_shift;
	std::vector<float> mono_data;
	std::vector<float> mono_state(mono_delay, 0.0);
	
	std::vector<float> left_data;
	std::vector<float> right_data;
	std::vector<float> stereo_data;
	
	std::vector<short int> audio;
	
	for (;; block_count++){
		
		readStdinBlockData(block_size, block_count, iq_block);
		if((std::cin.rdstate()) != 0){
			std::cerr << "End of input stream reached!" << std::endl;
			exit(1);
		}
		//std::cerr << "Read block " << block_count << std::endl;
		
		////////////////////////////////RF FRONT END////////////////////////////
		
		// separate iq_block into I and Q
		int j = 0;
		for (int i = 0; i < (int)iq_block.size(); i+=2){
			i_block[j] = iq_block[i];
			q_block[j] = iq_block[i+1];
			j++;
		};
		
		
		
		// resampling filter (2.4MS/s -> 240kS/s, 0 - 100kHz)
		resample(i_ds, state_i_lpf_100k, i_block, rf_coeff, 1, rf_decim);
		resample(q_ds, state_q_lpf_100k, q_block, rf_coeff, 1, rf_decim);
		
		// FM demodulation
		FMDemod(demod_data, prev_i, prev_q, i_ds, q_ds);
		
		////////////////////////////////RF FRONT END////////////////////////////
			
		monoProcessing(mono_data, audio_coeff, audio_state, demod_data, audio_decim, audio_interp, block_count);
				
		
		if (channels == 2){
			// MONO DELAY
			mono_shift.clear();
			mono_shift.reserve(mono_data.size());
			mono_shift.insert(mono_shift.end(), mono_state.begin(), mono_state.end());
			mono_shift.insert(mono_shift.end(), mono_data.begin(), mono_data.end() - mono_delay);
			
			mono_state.clear();
			mono_state.insert(mono_state.end(), mono_data.end() - mono_delay, mono_data.end());
			
			stereoProcessing(stereo_data, 
							carrier_coeff, carrier_state, channel_coeff, channel_state, audio_coeff, audio_state, 
							demod_data, audio_decim, audio_interp, bp_fs, 
							integrator, phaseEst, feedbackI, feedbackQ, ncoOut_state, trigOffset, block_count);
			
			// stereo combiner
			LRExtraction(left_data, right_data, mono_shift, stereo_data);

			int j = 0;
			audio.clear();
			audio.reserve(left_data.size() + right_data.size());
			audio.resize(left_data.size() + right_data.size(), 0.0);
			for (unsigned int k=0; k < left_data.size(); k++){
				// RIGHT
				if (std::isnan(right_data[k])) audio[j] = 0;
				// prepare a block of audio data to be redirected to stdout at once
				else audio[j] = static_cast<short int>(right_data[k] * 16384);
				// LEFT
				if (std::isnan(left_data[k])) audio[j+1] = 0;
				// prepare a block of audio data to be redirected to stdout at once
				else audio[j+1] = static_cast<short int>(left_data[k] * 16384);				
				j+=2;
			}
			
			fwrite(&audio[0], sizeof(short int), audio.size(), stdout);
		}
		else if (channels == 1){
			audio.clear();
			audio.reserve(mono_data.size());
			audio.resize(mono_data.size(), 0.0);
			for (unsigned int k=0; k < mono_data.size(); k++){
				if (std::isnan(mono_data[k])) audio[k] = 0;
				// prepare a block of audio data to be redirected to stdout at once
				else audio[k] = static_cast<short int>(mono_data[k] * 16384);
			}
			
			fwrite(&audio[0], sizeof(short int), audio.size(), stdout);
		}
	}
	
	std::cerr << "done" << std::endl;
	return 0;
}

// cat ../data/samples_u8.raw | ./project | aplay -c 1 -f S16_LE -r 48000
// rtl_sdr -f 87.5M -s 2.4M - | ./project | aplay -c 1 -f S16_LE -r 48000
// gnuplot -e 'set terminal png size 1024,768' ../data/example.gnuplot > ../data/example.png
