/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "../include/filter.h"
#include <cmath>

// function to compute the impulse response "h" based on the sinc function
void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h)
{
	// allocate memory for the impulse response
	h.clear(); h.resize(num_taps, 0.0);

	float normCutoff = Fc/(Fs/2);
	float inverse_taps = 1/(float)num_taps;

	for (int i = 0; i < num_taps; i++) {
		if (i == (num_taps-1) * 0.5) {
			h[i] = normCutoff;
		}
		else {
			float denominator = PI * normCutoff * (i - ((num_taps-1) * 0.5));
			float numerator = sin(denominator);

			h[i] = normCutoff * (numerator / denominator);
		}

		h[i] = h[i] * std::pow(sin(i * PI * inverse_taps), 2);
	}
}

// function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h" in blocks
void LPFilter(std::vector<float> &output, 
			const std::vector<float> &input, 
			const std::vector<float> &coeff, 
			std::vector<float> &state,
			int block_count)
{
	int nTaps = (int)coeff.size();
	int blockSize = (int)input.size();
	
	// allocate memory for the output (filtered) data
	output.clear(); output.resize(blockSize, 0.0);
		
	// concatenate state, input
	std::vector<float> signal;
	signal.insert(signal.end(), state.begin(), state.end());
	signal.insert(signal.end(), input.begin(), input.end());

	// discrete convolution
	for (int n = 0; n < blockSize; n++){
		for (int k = 0; k < nTaps; k++){
			//float prevOut = output[n];
			output[n] += coeff[k] * signal[n-k + nTaps-1];
			
			// print each sumproduct in calculating output
			//if (block_count == 0 && n == 0) {std::cerr << output[n] << " = " << prevOut << " + " << coeff[k] << " * " << signal[n-k + nTaps-1] << std::endl;}
		}
	}
	
	// state saving
	state.clear();
	state.resize(nTaps - 1);
	int indexState = 0;
	for (int c = blockSize - nTaps+1; c < blockSize; c++){
		state[indexState] = input[c];
		indexState++;
	}
}

void demodFM(const std::vector<float> &i_ds, const std::vector<float> &q_ds, std::vector<float> &demod, float p_i, float p_q) {
	for (unsigned int i=0; i < i_ds.size(); i++) {
		float currI = i_ds[i];
		float currQ = q_ds[i];

		float derivI = currI - p_i;
		float derivQ = currQ - p_q;

		if ((pow(currI,2) + pow(currQ,2)) == 0) {
			demod.push_back(0);
		}
		else {
			float numerator = currI*derivI - currQ*derivQ;
			float denominator = pow(currI,2) + pow(currQ,2);
			demod.push_back(numerator/denominator);
		}

		p_i = currI;
		p_q = currQ;
	}
}

void downSample(const std::vector<float> &original, std::vector<float> &downsampled, int decim) {

	for (unsigned int i=0; i < original.size(); i=i+decim) {
		downsampled.push_back(original[i]);
	}
}
