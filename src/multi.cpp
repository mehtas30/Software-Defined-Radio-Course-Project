#include <iostream>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <vector>

#include "../include/dy4.h"
#include "../include/filter.h"
#include "../include/fourier.h"
#include "../include/genfunc.h"
#include "../include/iofunc.h"
#include "../include/logfunc.h"
#include <chrono>

#define RANDOM_BOUND 1000
#define RANDOM_SEED 0
#define ELEM_SIZE 100
#define TOTAL_ELEMS 10
int generated_elems = 0;

#define QUEUE_ELEMS 3

void RF_thread(std::vector<float> iq_block,std::vector<float> i_ds, std::vector<float> q_ds, 	std::vector<float> i_block, std::vector<float> q_block,	std::vector<float> state_i_lpf_100k,std::vector<float> state_q_lpf_100k,std::vector<float> rf_coeff,std::vector<float> demod_data,	float prev_i,float prev_q, int block_count){
		// separate iq_block into I and Q
        int rf_decim = 10;
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
		
		if (block_count == 0){
			std::cerr << "demod_data=" << std::endl;
			for (int i = 0; i < 30; i++){
				std::cerr << demod_data[i] << "\n";
			}
		}

}
void monoProcessing(std::vector<float> &mono_data,
					const std::vector<float> &audio_coeff, std::vector<float> &audio_state,
					const std::vector<float> &demod_data,
					int &audio_decim, int &audio_interp) 
{
				
	// resampling filter (240kS/s -> 48kS/s, 0 - 16kHz)
	resample(mono_data, audio_state, demod_data, audio_coeff, audio_interp, audio_decim);

}


void stereoProcessing(std::vector<float> &left_data, std::vector<float> &right_data, 
					const std::vector<float> &mono_data, 
					const std::vector<float> &carrier_coeff, std::vector<float> &carrier_state,
					const std::vector<float> &channel_coeff, std::vector<float> &channel_state,
					const std::vector<float> &audio_coeff, std::vector<float> &audio_state,
					const std::vector<float> &demod_data,
					const int &audio_decim, const int &audio_interp, const int &if_fs, const int &block_count)
{
	// channel extraction
	std::vector<float> channel_data;
	
	// carrier recovery
	std::vector<float> carrier_data;
	
	// stereo processing
	std::vector<float> mixer_data;
	std::vector<float> stereo_data;

	

	// stereo channel extraction
	//resample(channel_data, channel_state, demod_data, channel_coeff, 1, 1);
	
	// stereo carrier recovery
	//resample(carrier_data, carrier_state, demod_data, carrier_coeff, 1, 1);
	PLL(carrier_data, 19000, if_fs, 2, 0, 0.01);
	
	// stereo processing
	mixer(mixer_data, channel_data, carrier_data);
	//resample(stereo_data, audio_state, mixer_data, audio_coeff, audio_interp, audio_decim);
	
	// at this point stereo_data contains L - R
	
	// stereo combiner
	LRExtraction(left_data, right_data, mono_data, stereo_data);

	// COMBINE LEFT AND RIGHT INTO ONE VECTOR -> left_right_data

}



void audio_thread(std::vector<float> &mono_data,
					const std::vector<float> &audio_coeff, std::vector<float> &audio_state,
					const std::vector<float> &demod_data,
					int &audio_decim, int &audio_interp,
					std::vector<float> &carrier_coeff, std::vector<float> &carrier_state, \
					std::vector<float> &channel_coeff, std::vector<float> &channel_state, \
				    const int &if_fs, const int &block_count, int &channels, std::vector<float> &mono_temp,std::vector<float> &mono_state,int &mono_delay,	std::vector<float> &left_data, \
	                std::vector<float> &right_data, \
	                std::vector<short int> &audio) {

		monoProcessing(mono_data, audio_coeff, audio_state, demod_data, audio_decim,audio_interp);
                // MONO DELAY
		mono_data.clear();
		mono_data.reserve(mono_temp.size());
		mono_data.insert(mono_data.end(), mono_state.begin(), mono_state.end());
		mono_data.insert(mono_data.end(), mono_temp.begin(), mono_temp.end() - mono_delay);
		
		mono_state.clear();
		mono_state.insert(mono_state.end(), mono_temp.end() - mono_delay, mono_temp.end());
		
		if (channels == 2){
			stereoProcessing(left_data, right_data, mono_data, carrier_coeff, carrier_state, channel_coeff, channel_state, audio_coeff, audio_state, demod_data, audio_decim, audio_interp, if_fs, block_count);

		}
		
		if (channels == 1){
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
		else if (channels == 2){
			int j = 0;
			audio.clear();
			audio.reserve(left_data.size() + right_data.size());
			audio.resize(left_data.size() + right_data.size(), 0.0);
			for (unsigned int k=0; k < left_data.size(); k++){
				//// RIGHT
				//if (std::isnan(right_data[k])) audio[j] = 0;
				//// prepare a block of audio data to be redirected to stdout at once
				//else audio[j] = static_cast<short int>(right_data[k] * 16384);
				// LEFT
				if (std::isnan(left_data[k])) audio[j+1] = 0;
				// prepare a block of audio data to be redirected to stdout at once
				else audio[j+1] = static_cast<short int>(left_data[k] * 16384);				
				j+=2;
			}
			
			fwrite(&audio[0], sizeof(short int), audio.size(), stdout);
		}
	}
	






void my_producer(std::queue<std::vector<int>> &my_queue,  \
	std::mutex& my_mutex, \
	std::condition_variable& my_cvar,std::vector<float> i_ds, \
std::vector<float>& q_ds, \ 
std::vector<float>& i_block, \ 
std::vector<float>& q_block, \ 
std::vector<float> &state_i_lpf_100k, \
std::vector<float> &state_q_lpf_100k,\
std::vector<float> &rf_coeff, \
std::vector<float> &demod_data,\
float &prev_i, \
float &prev_q,\ 
int &block_count, \ 
int &block_size, \
std::vector<float> &iq_block)
{           
	for (;; block_count++){
        readStdinBlockData(block_size, block_count, iq_block);
        if ((std::cin.rdstate())!=0){
            break;
        }

		std::vector<int> elem;
		RF_thread(iq_block,i_ds,q_ds,i_block,q_block,state_i_lpf_100k,state_q_lpf_100k,rf_coeff,demod_data,prev_i,prev_q,block_count);
		generated_elems++;

		std::unique_lock<std::mutex> my_lock(my_mutex);
		while (my_queue.size()>=QUEUE_ELEMS) {
			my_cvar.wait(my_lock);
		}
		my_queue.push(elem);
		my_cvar.notify_one();
		my_lock.unlock();
	}
}

void my_consumer(std::queue<std::vector<int>> &my_queue,  \
	std::mutex& my_mutex, \
	std::condition_variable& my_cvar, std::vector<float> &mono_data,
					const std::vector<float> &audio_coeff, std::vector<float> &audio_state,
					const std::vector<float> &demod_data,
					int& audio_decim, int& audio_interp,
					std::vector<float> &carrier_coeff, std::vector<float> &carrier_state, \
					std::vector<float> &channel_coeff, std::vector<float> &channel_state, \
				    const int &if_fs, const int &block_count, int &channels, std::vector<float> &mono_temp,std::vector<float> &mono_state,int &mono_delay,	std::vector<float> &left_data, \
	                std::vector<float> &right_data, \
	                std::vector<short int> &audio)
{
	while (1)
	{
		std::unique_lock<std::mutex> my_lock(my_mutex);
		while (my_queue.empty()) {
			my_cvar.wait(my_lock);
		}
		std::vector<int> elem = my_queue.front();
		my_queue.pop();
		my_cvar.notify_one();
		my_lock.unlock();

		audio_thread(mono_data,
					audio_coeff,audio_state,
					demod_data,
					audio_decim, audio_interp,
					carrier_coeff, carrier_state, \
					channel_coeff, channel_state, \
				    if_fs, block_count, channels, mono_temp, mono_state, mono_delay,left_data, \
	                right_data, \
	                audio);
		if ((generated_elems == TOTAL_ELEMS) && my_queue.empty())
			break;
	}
}

int main(int argc, char* argv[])
{
	srand(RANDOM_SEED);
	// binary files can be generated through the
	// Python models from the "../model/" sub-folder
	const std::string in_fname = "../data/fm_demod_10.bin";
	std::vector<float> bin_data;
	readBinData(in_fname, bin_data);

	int mode = 0;
	int channels = 2;

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
	int audio_taps = 51;
	int audio_decim = 5;
	int audio_interp = 1;

	int block_size = 128 * rf_decim * audio_decim * 2;
	int block_count = 0;
	
	int mono_delay = (audio_taps - 1) * 0.5;

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
	std::vector<float> channel_state(audio_taps-1, 0.0);
	// coefficients for stereo channel extraction, 22kHz to 54kHz
	std::vector<float> channel_coeff;
	impulseResponseBPF(channel_coeff, if_fs, 22000.0, 54000.0, audio_taps);
	
	// stereo carrier recovery 
	std::vector<float> carrier_state(audio_taps-1, 0.0);
	// coefficients for stereo carrier recovery, 18.5kHz to 19.5kHz
	std::vector<float> carrier_coeff;
	impulseResponseBPF(carrier_coeff, if_fs, 18500, 19500, audio_taps);
	
	// state saving for extracting audio (Fc = 16kHz)
	std::vector<float> audio_state(audio_taps-1, 0.0);
	// coefficients for IF -> audio LPF, Fc = 16kHz
	std::vector<float> audio_coeff;
	impulseResponseLPF(audio_coeff, if_fs, audio_fc, audio_taps, audio_interp);
	
	std::vector<float> mono_temp;
	std::vector<float> mono_data;
	std::vector<float> mono_state(mono_delay, 0.0);
	
	std::vector<float> left_data;
	std::vector<float> right_data;
	std::vector<short int> audio;
	
	// for (;; block_count++){
		
	// 	readStdinBlockData(block_size, block_count, iq_block);
	// 	if((std::cin.rdstate()) != 0){
	// 		std::cerr << "End of input stream reached!" << std::endl;
	// 		exit(1);
	// 	}
    // }
		//std::cerr << "Read block " << block_count << std::endl;
		
	
	std::queue<std::vector<int>> my_queue;
	std::mutex my_mutex;
	std::condition_variable my_cvar;

	std::thread RF_producer = std::thread(my_producer, std::ref(my_queue), \
		std::ref(my_mutex), std::ref(my_cvar), i_ds, 
         q_ds, \ 
        i_block, \ 
        q_block, \ 
        state_i_lpf_100k, \
        state_q_lpf_100k,\
        rf_coeff, \
        demod_data,\
        prev_i, \
        prev_q,\ 
        block_count, \ 
        block_size, \
        iq_block);

	std::thread audio_consumer = std::thread(my_consumer, std::ref(my_queue), \
		std::ref(my_mutex), mono_data,
					audio_coeff,audio_state,
					demod_data,
					audio_decim, audio_interp,
					carrier_coeff, carrier_state, \
					channel_coeff, channel_state, \
				    if_fs, block_count, channels, mono_temp, mono_state, mono_delay,left_data, \
	                right_data, \
	                audio);

	RF_producer.join();
	audio_consumer.join();

	return 0;
}