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

void RF_thread(std::vector<float> &mono_data,
					const std::vector<float> &audio_coeff, std::vector<float> &audio_state,
					const std::vector<float> &demod_data,
					int audio_decim, int audio_interp,std::queue<std::vector<int>> &my_queue,  \
	                std::mutex& my_mutex, \
	                std::condition_variable& my_cvar){

}



void audio_thread(std::vector<float> &mono_data,
					const std::vector<float> &audio_coeff, std::vector<float> &audio_state,
					const std::vector<float> &demod_data,
					int audio_decim, int audio_interp,std::queue<std::vector<int>> &my_queue,  \
	                std::mutex& my_mutex, \
	                std::condition_variable& my_cvar) {

				

}


void my_producer(std::queue<std::vector<int>> &my_queue,  \
	std::mutex& my_mutex, \
	std::condition_variable& my_cvar)
{
	while (generated_elems < TOTAL_ELEMS)
	{
		std::vector<int> elem;
		RF_thread();
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

void my_consume(std::queue<std::vector<int>> &my_queue,  \
	std::mutex& my_mutex, \
	std::condition_variable& my_cvar)
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

		audio_thread();
		if ((generated_elems == TOTAL_ELEMS) && my_queue.empty())
			break;
	}
}

int main()
{
	srand(RANDOM_SEED);

	std::queue<std::vector<int>> my_queue;
	std::mutex my_mutex;
	std::condition_variable my_cvar;

	std::thread RF_producer = std::thread(RF_thread, std::ref(my_queue), \
		std::ref(my_mutex), std::ref(my_cvar));

	std::thread audio_consumer = std::thread(audio_thread, std::ref(my_queue), \
		std::ref(my_mutex), std::ref(my_cvar));

	RF_producer.join();
	audio_consumer.join();

	return 0;
}