/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Copyright by Nicola Nicolici
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>
#include "../include/dy4.h"
#include "../include/filter.h"
#include "../include/iofunc.h"

#define QUEUE_CAPACITY 3

void rf_thread(std::queue<std::vector<float>> &queue, std::queue<std::vector<float>> &rds_queue,
			   std::mutex &mutex, std::mutex &rds_mutex, 
			   std::condition_variable &cvar, std::condition_variable &rds_cvar,
			   const int &rf_fs, const int &rf_fc, const int &rf_taps, const int &rf_decim, 
			   const int &block_size, int &block_count)
{
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
	
	while(1){
		
		readStdinBlockData(block_size, block_count, iq_block);
		if((std::cin.rdstate()) != 0){
			std::cerr << "End of input stream reached!" << std::endl;
			exit(1);
		}
		
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

		std::unique_lock<std::mutex> lock(mutex);
		//std::unique_lock<std::mutex> rds_lock(rds_mutex);
		while(queue.size() >= QUEUE_CAPACITY){
			cvar.wait(lock);
			//rds_cvar.wait(rds_lock);
		}
		
		queue.push(demod_data);
		cvar.notify_one();
		lock.unlock();
		//rds_queue.push(demod_data);
		//rds_cvar.notify_one();
		//rds_lock.unlock();
	}
}

void audio_thread(std::queue<std::vector<float>> &queue,
				  std::mutex &mutex, std::condition_variable &cvar,
				  const int &bp_fs, const int &bp_taps, const int &if_fs, const int &mono_delay,
				  const int &audio_fc, const int &audio_taps, const int &audio_decim, const int &audio_interp)
{
	// stereo channel extraction
	std::vector<float> channel_data;
	std::vector<float> channel_state(bp_taps-1, 0.0);
	// coefficients for stereo channel extraction, 22kHz to 54kHz
	std::vector<float> channel_coeff;
	impulseResponseBPF(channel_coeff, bp_fs, 22000.0, 54000.0, bp_taps);
	
	// stereo carrier recovery 
	std::vector<float> carrier_data;
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
	
	// stereo processing
	std::vector<float> mixer_data;
	
	std::vector<float> left_data;
	std::vector<float> right_data;
	std::vector<float> stereo_data;
	
	std::vector<short int> audio;
	
	while(1){
		std::unique_lock<std::mutex> lock(mutex);
		while(queue.empty()){
			cvar.wait(lock);
		}
		
		std::vector<float> demod_data = queue.front();
		queue.pop();
		cvar.notify_one();
		lock.unlock();
		
		//////////////////////////////MONO//////////////////////////////
		
		// resampling filter (240kS/s -> 48kS/s, 0 - 16kHz)
		resample(mono_data, audio_state, demod_data, audio_coeff, audio_interp, audio_decim);
		
		//////////////////////////////MONO//////////////////////////////
		
		/////////////////////////////STEREO/////////////////////////////
		
		// mono delay
		mono_shift.clear();
		mono_shift.reserve(mono_data.size());
		mono_shift.insert(mono_shift.end(), mono_state.begin(), mono_state.end());
		mono_shift.insert(mono_shift.end(), mono_data.begin(), mono_data.end() - mono_delay);
		
		mono_state.clear();
		mono_state.insert(mono_state.end(), mono_data.end() - mono_delay, mono_data.end());
		
		// stereo channel extraction
		resample(channel_data, channel_state, demod_data, channel_coeff, 1, 1);

		// stereo carrier recovery
		resample(carrier_data, carrier_state, demod_data, carrier_coeff, 1, 1);
		PLL(carrier_data, 19000, if_fs, 2, 0, 0.01, integrator, phaseEst, feedbackI, feedbackQ, ncoOut_state, trigOffset);
		
		// stereo processing
		mixer(mixer_data, channel_data, carrier_data);
		
		// 16kHz filter
		resample(stereo_data, audio_state, mixer_data, audio_coeff, audio_interp, audio_decim);
		
		// stereo combiner
		LRExtraction(left_data, right_data, mono_shift, stereo_data);
		
		/////////////////////////////STEREO/////////////////////////////
		
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
}


void rds_thread(std::queue<std::vector<float>> &queue,
				std::mutex &mutex, std::condition_variable &cvar,
				const int &bp_fs, const int &bp_taps, const int &rds_fs, const int &channel_delay,
				const int &rds_fc, const int &rds_taps, const int &rds_decim, const int &rds_interp)
{
	// channel extraction
	std::vector<float> channel_shift;
	std::vector<float> channel_shift_state(channel_delay, 0.0);
	std::vector<float> channel_data;
	std::vector<float> channel_state(bp_taps-1, 0.0);
	std::vector<float> extract_coeff;
	impulseResponseBPF(extract_coeff, bp_fs, 54000, 60000, bp_taps);
	
	// carrier recovery
	std::vector<float> channel_squared;
	std::vector<float> carrier_data;
	std::vector<float> carrier_state(bp_taps-1, 0.0);
	std::vector<float> carrier_coeff;
	impulseResponseBPF(carrier_coeff, bp_fs, 113500, 114500, bp_taps);
	// PLL state saving
	float integrator = 0.0;
	float phaseEst = 0.0;
	float feedbackI = 1.0;
	float feedbackQ = 0.0;
	float trigOffset = 0.0;
	float ncoOut_state = 1.0;
	
	std::vector<float> mixer_data;
	
	// 3kHz LPF
	std::vector<float> demod_coeff;
	impulseResponseLPF(demod_coeff, bp_fs, rds_fc, rds_taps, 1);
	
	while(1){
		std::unique_lock<std::mutex> lock(mutex);
		while(queue.empty()){
			cvar.wait(lock);
		}
		
		std::vector<float> demod_data = queue.front();
		queue.pop();
		cvar.notify_one();
		lock.unlock();
		
		// channel extraction
		resample(channel_data, channel_state, demod_data, extract_coeff, 1, 1);
		
		// carrier recovery
		channel_squared.clear(); channel_squared.reserve(channel_data.size());
		channel_squared.resize(channel_data.size(), 0.0);
		for (int i = 0; i < channel_data.size(); i++){
			channel_squared[i] = channel_data[i] * channel_data[i];
		}
		
		// 113.5kHz - 114.5kHz BPF
		resample(carrier_data, carrier_state, channel_squared, carrier_coeff, 1, 1);
		
		PLL(carrier_data, 114000, bp_fs, 0.5, 0, 0.01, integrator, phaseEst, feedbackI, feedbackQ, ncoOut_state, trigOffset);
		
		// channel delay
		channel_shift.clear();
		channel_shift.reserve(channel_data.size());
		channel_shift.insert(channel_shift.end(), channel_shift_state.begin(), channel_shift_state.end());
		channel_shift.insert(channel_shift.end(), channel_data.begin(), channel_data.end() - channel_delay);
		
		channel_shift_state.clear();
		channel_shift_state.insert(channel_shift_state.end(), channel_data.end() - channel_delay, channel_data.end());
		
		// RDS demodulation mixer
		mixer(mixer_data, carrier_data, channel_shift);
	}
}

int main(int argc, char* argv[])
{
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
	
	const int rf_fc = 100000;
	const int audio_fc = 16000;
	const int rf_taps = 51;
	const int bp_taps = 51;
	const int mono_delay = 5;
	// const int rds_fc = 3000;
	// const int rds_delay = 5;

	int rf_fs = 2400000;
	int rf_decim = 10;
	
	int if_fs = 240000;
	int bp_fs = 240000;
	int rds_fs = 240000;

	int audio_taps = 51;
	int audio_decim = 5;
	int audio_interp = 1;

	//int rds_taps = 51;
	//int rds_decim = 384;
	//int rds_interp = 19;

	if (mode == 0){
		rf_fs = 2400000; rf_decim = 10;
		if_fs = 240000; bp_fs = 240000;
		audio_decim = 5;
		//rds_taps *= rds_interp;
		//rds_fs *= rds_interp;
	}
	else if (mode == 1){ 
		rf_fs = 1152000; rf_decim = 4;
		if_fs = 288000; bp_fs = 288000; //rds_fs = 288000;
		audio_decim = 6;
		//rds_decim = 2304;
		//rds_interp = 95;
		//rds_taps *= rds_interp;
		//rds_fs *= rds_interp;
	}
	else if (mode == 2){ 
		rf_fs = 2400000; rf_decim = 10;
		if_fs = 240000; bp_fs = 240000;
		audio_decim = 800; audio_interp = 147;
		audio_taps *= audio_interp;
		if_fs *= audio_interp;
		//rds_taps *= rds_interp;
		//rds_fs *= rds_interp;
	}
	else if (mode == 3){ 
		rf_fs = 2304000; rf_decim = 9;
		if_fs = 256000; bp_fs = 256000; rds_fs = 256000;
		audio_decim = 2560; audio_interp = 441;
		audio_taps *= audio_interp;
		if_fs *= audio_interp;
		//rds_decim = 2048;
		//rds_interp = 95;
		//rds_taps *= rds_interp;
		//rds_fs *= rds_interp;
	}
	
	int block_size = 256 * rf_decim * audio_decim;
	int block_count = 0;

	std::queue<std::vector<float>> queue;
	std::queue<std::vector<float>> rds_queue;
	std::mutex mutex;
	std::mutex rds_mutex;
	std::condition_variable cvar;
	std::condition_variable rds_cvar;
	
	std::thread t_rf = std::thread(rf_thread, std::ref(queue), std::ref(rds_queue), std::ref(mutex), std::ref(rds_mutex), std::ref(cvar), std::ref(rds_cvar),
								   std::ref(rf_fs), std::ref(rf_fc), std::ref(rf_taps), std::ref(rf_decim), 
								   std::ref(block_size), std::ref(block_count));
	std::thread t_audio = std::thread(audio_thread, std::ref(queue), std::ref(mutex), std::ref(cvar),
									  std::ref(bp_fs), std::ref(bp_taps), std::ref(if_fs), std::ref(mono_delay),
									  std::ref(audio_fc), std::ref(audio_taps), std::ref(audio_decim), std::ref(audio_interp));
	//std::thread t_rds = std::thread(rds_thread, std::ref(rds_queue), std::ref(rds_mutex), std::ref(rds_cvar),
									//std::ref(bp_fs), std::ref(bp_taps), std::ref(rds_fs), std::ref(rds_delay),
									//std::ref(rds_fc), std::ref(rds_taps), std::ref(rds_decim), std::ref(rds_interp));
	
	t_rf.join();
	t_audio.join();
	// t_rds.join();

	std::cerr << "done" << std::endl;
	return 0;
}

// cat ../data/samples_u8.raw | ./project | aplay -c 1 -f S16_LE -r 48000
// rtl_sdr -f 102.9M -s 2.4M - | ./project 0 2 | aplay -c 2 -f S16_LE -r 48000
// gnuplot -e 'set terminal png size 1024,768' ../data/example.gnuplot > ../data/example.png
