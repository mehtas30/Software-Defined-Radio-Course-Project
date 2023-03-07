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

	// the rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab
	float normCutoff = (Fc*2)/Fs;

	for (int i = 0; i < num_taps; i++) {
		if (i == (num_taps-1)/2) {
			h[i] = 0;
		}
		else {
			float numerator = sin(PI * normCutoff * (i - (num_taps - 1)/2));
			float denominator = PI * normCutoff * (i - (num_taps-1)/2);

			h[i] = normCutoff * (numerator/denominator);
		}

		h[i] = h[i] * pow((sin((i * PI)/num_taps)), 2);
	}
}

// function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h"
void convolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h)
{
	// allocate memory for the output (filtered) data
	y.clear(); y.resize(x.size(), 0.0);

	// the rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab
	for (int n = 0; n < (int)(x.size()); n++){
		for (int k = 0; k < (int)(h.size()); k++){
			if (n-k >= 0 && n-k < (int)(x.size())){
				y[n] += h[k] * x[n-k];
			}
		}
	}
}
