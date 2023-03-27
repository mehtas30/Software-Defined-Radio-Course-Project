/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "../include/filter.h"
#include <cmath>
#include <chrono>

// function to compute the impulse response "h" based on the sinc function
void impulseResponseLPF(std::vector<float> &h, const float Fs, const float Fc, const int num_taps, const int gain)
{
	// allocate memory for the impulse response
	h.clear(); h.reserve(num_taps); h.resize(num_taps, 0.0);

	float norm_fc = Fc / (Fs / 2);
	float inverse_taps = 1 / (float)num_taps;

	for (int i = 0; i < num_taps; i++) {
		if (i == (num_taps - 1) * 0.5) {
			h[i] = norm_fc;
		}
		else {
			float denominator = PI * norm_fc * (i - ((num_taps - 1) * 0.5));
			float numerator = sin(denominator);

			h[i] = norm_fc * (numerator / denominator);
		}

		h[i] *= std::pow(sin(i * PI * inverse_taps), 2);

		if (gain != 1){ h[i] *= gain; }
	}
}

// function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h" in blocks
void LPFilter(std::vector<float> &output,
			std::vector<float> &state,
			const std::vector<float> &input,
			const std::vector<float> &coeff)
{
	int taps = (int)coeff.size();
	int block_size = (int)input.size();
	int state_size = (int)state.size();

	//std::cerr << "taps=" << taps << ", block size=" << block_size << std::endl;
	// allocate memory for the output (filtered) data
	output.clear(); 
	output.reserve(block_size);
	
	//std::cerr << "reserved size = " << output.size() << std::endl;
		
	// concatenate state, input
	std::vector<float> signal;
	signal.insert(signal.end(), state.begin(), state.end());
	signal.insert(signal.end(), input.begin(), input.end());

	// discrete convolution
	for (int n = 0; n < block_size; n++){
		float sum_product = 0.0;
		for (int k = 0; k < taps; k++){
			int nk = n-k;
			if ((nk >= 0) && (nk < block_size)){
				sum_product += coeff[k] * input[nk];
			}
			else {
				sum_product += coeff[k] * state[state_size + nk];
			}
		}
		output.push_back(sum_product);
	}
	
	//std::cerr << "filled size = " << output.size() << std::endl;
	// state saving
	state.clear();
	state.resize(taps - 1);
	int indexState = 0;
	for (int c = block_size - taps+1; c < block_size; c++){
		state[indexState] = input[c];
		indexState++;
	}
}

void resample(std::vector<float> &output,
			std::vector<float> &state,
			const std::vector<float> &input,
			const std::vector<float> &coeff,
			const int up_factor,
			const int down_factor)
{
	int taps = (int)coeff.size();
	int block_size = (int)input.size();
	int state_size = (int)state.size();

	//std::cerr << "taps=" << taps << ", block size=" << block_size << std::endl;
	// allocate memory for the output (filtered) data
	output.clear(); 
	output.reserve(block_size * up_factor / down_factor);
	output.resize(block_size * up_factor / down_factor, 0.0);
	
	for (int n = 0; n < output.size(); n++){
		int phase = (n * down_factor) % up_factor;
		
		for (int k = phase; k < taps; k += up_factor){
			int j = int(((down_factor * n) - k) / up_factor);
			
			if ((j >= 0) && (j < block_size)){
				output[n] += coeff[k] * input[j];
			} 
			else{
				output[n] += coeff[k] * state[state_size + j];
			}
		}
	}
	
	// state saving
	state.clear();
	state.resize(taps - 1);
	int indexState = 0;
	for (int c = block_size - taps + 1; c < block_size; c++){
		state[indexState] = input[c];
		indexState++;
	}
}


void FMDemod(std::vector<float> &fm_demod, float &prev_i, float &prev_q, const std::vector<float> &i_ds, const std::vector<float> &q_ds) {
	
	fm_demod.clear(); fm_demod.reserve(i_ds.size());
	
	for (int i = 0; i < (int)i_ds.size(); i++) {
		float curr_i = i_ds[i];
		float curr_q = q_ds[i];

		float deriv_i = curr_i - prev_i;
		float deriv_q = curr_q - prev_q;
		
		// i^2 + q^2
		float denominator = std::pow(curr_i, 2) + std::pow(curr_q, 2);

		if (denominator != 0) {
			// i * dq/dt - q * di/dt
			float numerator = (curr_i * deriv_q) - (curr_q * deriv_i);
			fm_demod.push_back(numerator / denominator);
		}

		// state saving
		prev_i = curr_i;
		prev_q = curr_q;
	}
}

void downsample(std::vector<float> &downsampled, const std::vector<float> &data, const int down_factor) {

	downsampled.clear(); downsampled.reserve(data.size());
	
	for (int i = 0; i < data.size(); i += down_factor) {
		//std::cerr << "c0" << std::endl;
		downsampled.push_back(data[i]);
		//std::cerr << "c1" << std::endl;
	}
}

void upsample(std::vector<float> &upsampled, const std::vector<float> &data, const int up_factor){
	upsampled.clear();
	if (up_factor == 1) {
		upsampled.insert(upsampled.end(), data.begin(), data.end());
		return;
	}
	
	upsampled.reserve(data.size() * up_factor);
	upsampled.resize(data.size() * up_factor, 0.0);

	for (int i = 0; i < (int)data.size(); i++){
		upsampled[i * up_factor] = data[i];
	}
}

void bandpassfilter(std::vector<float> &output, float fb, float fe, float fs, int n_taps)
{
	float normCent = ((fe + fb) / 2) / (fs / 2);
	float normPass = (fe - fb) / (fs / 2);
	output.resize(n_taps - 1, 0.0);
	float value;
	for (int i = 0; i < (n_taps)-1; i++)
	{
		if (i == ((n_taps - 1) / 2))
		{
			output[i] = normPass;
		}
		else
		{
			value = PI * (normPass / 2) * (i - ((n_taps - 1) / 2));
			output[i] = normPass * sin(value) / value;
		}
		output[i] = output[i] * cos(i * PI * normCent);
		output[i] = output[i] * std::pow(sin(i * PI / n_taps), 2);
	}
}

void fmPLL(std::vector<float> &ncoOut, float freq, float Fs, float nocoScale, float phaseAdjust, float normBandwidth)
{
	float Cp = 2.666;
	float Ci = 3.555;
	float Kp = normBandwidth * Cp;
	float Ki = normBandwidth * normBandwidth * Ci;
	std::vector<float> pllin;
	pllin = ncoOut;
	ncoOut.clear();
	ncoOut.resize(pllin.size());
	float integrator = 0.0;
	float phaseEst = 0.0;
	float feedbackI = 0.0;
	float feedbackQ = 0.0;
	ncoOut[0] = 1.0;
	float trigOffset = 0.0;
	float errorI;
	float errorQ;
	float errorD;
	float trigArg;
	for (int i = 0; i < pllin.size(); i++)
	{
		errorI = pllin[i] * feedbackI;
		errorQ = pllin[i] * (-feedbackQ);
		errorD = atan2(errorQ, errorI);
		integrator = integrator + (Ki * errorD);
		phaseEst = phaseEst + (Kp * errorD) + integrator;
		trigOffset += 1;
		trigArg = 2 * PI * (freq / Fs) * (trigOffset) + phaseEst;
		feedbackI = cos(trigArg);
		feedbackQ = sin(trigArg);
		ncoOut[i] = cos(trigArg * nocoScale + phaseAdjust);
	}
}

void mixer(std::vector<float> &out, std::vector<float> &arr1, std::vector<float> &arr2)
{
	out.clear();
	for (int i = 0; i < arr1.size(); i++)
	{
		out.push_back(arr1[i] * arr2[i]);
	}
}

void lrExtraction(std::vector<float> &left, std::vector<float> &right, std::vector<float> monoData, std::vector<float> stereoData)
{
	left.clear();
	right.clear();
	left.resize(monoData.size());
	right.resize(monoData.size());
	for (int i = 0; i < monoData.size(); i++)
	{
		left[i] = (monoData[i] + stereoData[i]) / 2;
		right[i] = (monoData[i] - stereoData[i]) / 2;
	}
}
