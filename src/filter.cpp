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

	float normCutoff = (Fc*2)/Fs;
	float inverse_taps = 1/(float)num_taps;

	for (int i = 0; i < num_taps; i++) {
		if (i == (num_taps-1) * 0.5) {
			h[i] = 0;
		}
		else {
			float numerator = sin(PI * normCutoff * (i - ((num_taps-1) * 0.5)));
			float denominator = PI * normCutoff * (i - ((num_taps-1) * 0.5));

			h[i] = normCutoff * (numerator / denominator);
		}

		h[i] = h[i] * std::pow(sin(i * PI * inverse_taps), 2);
	}
}

// function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h" in blocks
void LPFilter(std::vector<float> &y, 
			const std::vector<float> &x, 
			const std::vector<float> &h, 
			std::vector<float> &state)
{
	// allocate memory for the output (filtered) data
	y.clear(); y.resize(x.size(), 0.0);

	// discrete convolution
	for (int n = 0; n < (int)(x.size()); n++){
		for (int k = 0; k < (int)(h.size()); k++){
			if (n-k >= 0){
				y[n] += h[k] * x[n-k];
			}
			// negative n-k correspond to right end of prev block
			else{ y[n] += h[k] * state[n-k]; }
		}
		y[n] = y[n] * 10;
	}
	// state saving
	std::copy(x.end()-state.size(), x.end(), state.begin());
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
