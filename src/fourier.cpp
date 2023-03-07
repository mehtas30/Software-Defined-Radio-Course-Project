/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

// source code for Fourier-family of functions
#include "dy4.h"
#include "fourier.h"

// just DFT function (no FFT yet)
void DFT(const std::vector<float> &x, std::vector<std::complex<float>> &Xf) {
	Xf.clear(); Xf.resize(x.size(), static_cast<std::complex<float>>(0));
	for (int m = 0; m < (int)Xf.size(); m++) {
		for (int k = 0; k < (int)x.size(); k++) {
				std::complex<float> expval(0, -2*PI*(k*m) / x.size());
				Xf[m] += x[k] * std::exp(expval);
		}
	}
}

// function to compute the magnitude values in a complex vector
void computeVectorMagnitude(const std::vector<std::complex<float>> &Xf, std::vector<float> &Xmag)
{
	// only the positive frequencies
	Xmag.clear(); Xmag.resize(Xf.size(), static_cast<float>(0));
	for (int i = 0; i < (int)Xf.size(); i++) {
		Xmag[i] = std::abs(Xf[i])/Xf.size();
	}
}

// function to estimate the PSD
void estimatePSD(std::vector<float> &freq,
				std::vector<float> &psd_est,
				const std::vector<float> &samples,
				const int freq_bins,
				const float Fs)
{
	// frequency increment

	float df = Fs / freq_bins;
	int freq_size = (int)(freq_bins / 2);
	int half_freq_bins = (int)(freq_bins / 2);
	psd_est.clear();
	psd_est.resize(half_freq_bins, 0.0);
	// create frequency vector to be used on the X axis
	freq.clear(); freq.resize(freq_size, 0.0);
	for (int i = 0; i < freq_size; i++){
		freq[i] = i*df;
	}

	// Hann window used to smooth the discrete data in order
	// to reduce spectral leakage after Fourier Transform
	std::vector<float> hann(freq_bins);
	for (int i = 0; i < freq_bins; i++){
		hann[i] = std::pow(std::sin(i * PI / freq_bins), 2);
	}

	// empty vector where PSD for each segment is computed
	std::vector<float> psd_list;

	/* samples should be a multiple of frequency bins, so
	number of segments used for estimation is an integer */
	int num_segments = (int)(samples.size() / freq_bins);

	std::vector<float> samples_slice(freq_bins, 0.0);
	std::vector<float> windowed_samples(freq_bins, 0.0);
	samples_slice.clear();

	// iterate through segments
	for (int k = 0; k < num_segments; k++){
		samples_slice.clear();
		windowed_samples.clear();
		windowed_samples.resize(freq_bins,0.0);
		auto start = samples.begin() + (k*freq_bins);
		auto end = samples.begin() + ((k+1)*freq_bins);
		copy(start, end, samples_slice.begin());

		// apply Hann window (point-wise multiplication)
		for (int i = 0; i < freq_bins; i++){
			windowed_samples[i] = samples_slice[i] * hann[i];
			
		}

		std::vector<std::complex<float>> Xf;
		// compute Fourier transform
		DFT(windowed_samples, Xf);

		std::vector<std::complex<float>> Xf_positive(half_freq_bins,static_cast<std::complex<float>>(0));
		// keep positive half of spectrum since input is real
		copy(Xf.begin(), Xf.begin() + half_freq_bins, Xf_positive.begin());
		std::vector<float> psd_seg(half_freq_bins, 0.0);

		for (int i = 0; i < (int)(psd_seg.size()); i++){
			// compute signal power and double for negative frequency half
			psd_seg[i] = (4 / (Fs * freq_bins)) * std::abs((std::pow(Xf_positive[i], 2)));
			
			// translate to dB
			psd_seg[i] = 10 * std::log10(psd_seg[i]);
			
		}
		// append
		psd_list.insert(psd_list.end(), psd_seg.begin(), psd_seg.end());
	}

	// iterate through positive frequency bins and average them
	for (int k = 0; k < half_freq_bins; k++){
		for (int l = 0; l < num_segments; l++){
			psd_est[k] += psd_list[k + (l * half_freq_bins)];
		}

		// compute estimate for each bin
		psd_est[k] = psd_est[k] / num_segments;
	}
}

// added IDFT
void IDFT(const std::vector<std::complex<float>> &Xf, std::vector<std::complex<float>> &x) {
	x.clear(); x.resize(Xf.size(), static_cast<std::complex<float>>(0));
	for (int k = 0; k < (int)x.size(); k++) {
		for (int m = 0; m < (int)x.size(); m++) {
			std::complex<float> expval(0, 2*PI*(k*m) / Xf.size());
			x[k] += Xf[m] * std::exp(expval);
		}
		x[k] /= Xf.size();
	}
}

// added FFT

// this is a support function used for bit reversal
unsigned int swap_bits(unsigned int x, unsigned char i, unsigned char j) {

	unsigned char bit_i = (x >> i) & 0x1L;
	unsigned char bit_j = (x >> j) & 0x1L;

	unsigned int val = x;
	val = bit_i ? (val | (0x1L << j)) : (val & ~(0x1L << j));
	val = bit_j ? (val | (0x1L << i)) : (val & ~(0x1L << i));

	return val;
}

// bit reversal is needed for the non-recursive (forward traversal) version of FFT
unsigned int bit_reversal(unsigned int x, unsigned char bit_size) {

	unsigned int val = x;

	for (int i=0; i < int(bit_size/2); i++)
		val = swap_bits(val, i, bit_size-1-i);

	return val;
}

// precomputation of twiddle factors (needed for faster versions of FFT)
void compute_twiddles(std::vector<std::complex<float>> &twiddles) {
	for (int k=0; k<(int)twiddles.size(); k++) {
			std::complex<float> expval(0.0, -2*PI*float(k)/ NFFT);
			twiddles[k] = std::exp(expval);
	}
}

// recursive version of FFT
// based on first principles from the mathematical derivation
void FFT_recursive(const std::vector<std::complex<float>> &x, \
	std::vector<std::complex<float>> &Xf) {

	if (x.size() > 1) {
		// declare vectors and allocate space for the even and odd halves
		std::vector<std::complex<float>> xe(int(x.size()/2)), xo(int(x.size()/2));
		std::vector<std::complex<float>> Xfe(int(x.size()/2)), Xfo(int(x.size()/2));

		// split into even and odd halves
		for (int k=0; k<(int)x.size(); k++)
			if ((k%2) == 0) xe[k/2] = x[k];
			else xo[k/2] = x[k];

		// call recursively FFT of half size for even and odd halves respectively
		FFT_recursive(xe, Xfe);
		FFT_recursive(xo, Xfo);

		// merge the results from the odd/even FFTs (each of half the size)
		for (int k=0; k<(int)xe.size(); k++) {
				std::complex<float> expval(0.0, -2*PI*float(k)/ x.size());
				std::complex<float> twiddle = std::exp(expval);
				Xf[k]           = Xfe[k] + twiddle * Xfo[k];
				Xf[k+xe.size()] = Xfe[k] - twiddle * Xfo[k];
		}
	} else {
		// end of recursion - copy time domain samples to frequency bins (default values)
		Xf[0] = x[0];
	}
}

// improved version of FFT
// the speed-up is mainly from avoiding the recomputation of twiddle factors
// the math library functions do not need to be recalled at each level of recursion
void FFT_improved(const std::vector<std::complex<float>> &x, \
	std::vector<std::complex<float>> &Xf, \
	const std::vector<std::complex<float>> &twiddles, \
	const unsigned char recursion_level) {

	if (x.size() > 1) {
		int half_size = int(x.size()/2);
		std::vector<std::complex<float>> xe(half_size), xo(half_size);
		std::vector<std::complex<float>> Xfe(half_size), Xfo(half_size);

		for (int k=0; k<half_size; k++) {
			xe[k] = x[k*2];
			xo[k] = x[k*2+1];
		}

		FFT_improved(xe, Xfe, twiddles, recursion_level+1);
		FFT_improved(xo, Xfo, twiddles, recursion_level+1);

		for (int k=0; k<half_size; k++) {
				// the precomputed twiddle factors are reused at each recursion level
				Xf[k]           = Xfe[k] + twiddles[k*(1<<(recursion_level-1))] * Xfo[k];
				Xf[k+half_size] = Xfe[k] - twiddles[k*(1<<(recursion_level-1))] * Xfo[k];
		}
	} else {
		Xf[0] = x[0];
	}
}

// optimized version of FFT is based on non-recursive forward traversal
// speed-up comes (in-part) from avoiding loading/unloading the stack with recursive calls
// however, in addition to pre-computation of twiddle factors, the central reason for
// speed-up is the in-place computation of the intermediate vectors at each FFT level
// the core idea is that ONLY two input values are used to compute two output levels
void FFT_optimized(const std::vector<std::complex<float>> &x, \
	std::vector<std::complex<float>> &Xf, \
	const std::vector<std::complex<float>> &twiddles) {

	unsigned char no_levels = (unsigned char)std::log2((float)x.size());
	for (int i=0; i<(int)x.size(); i++) {
		// time-domain values indexed in a bit-reversed manner before
		// copied to the leftmost (first) level from where forward traversal starts
		Xf[i] = x[bit_reversal(i, no_levels)];
	}

	unsigned int step_size = 1;

	std::complex<float> tmp;
	for (unsigned char l=0; l<no_levels; l++) {
		for (unsigned int p=0; p<x.size(); p+=2*step_size) {
			for (unsigned int k=p; k<p+step_size; k++) {
				// in-place computation - only two output values depend on two input values
				// indexing of twiddle factors exploits the fact that
				// index k for FFT of size n is the same as index 2k for FFT of size 2n
				tmp             = Xf[k] + twiddles[(k-p)*(1<<(no_levels-1-l))] * Xf[k+step_size];
				Xf[k+step_size] = Xf[k] - twiddles[(k-p)*(1<<(no_levels-1-l))] * Xf[k+step_size];
				Xf[k]           = tmp;
			}
		}
		step_size *= 2;
	}
}
