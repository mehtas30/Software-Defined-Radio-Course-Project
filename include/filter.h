/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_FILTER_H
#define DY4_FILTER_H

// add headers as needed
#include <iostream>
#include <vector>

// declaration of a function prototypes
void impulseResponseLPF(std::vector<float> &, const float, const float, const int, const int);

void LPFilter(std::vector<float> &, std::vector<float> &, const std::vector<float> &, const std::vector<float> &);

void resample(std::vector<float> &output, std::vector<float> &state, const std::vector<float> &input, const std::vector<float> &coeff, const int up_factor, const int down_factor);

void FMDemod(std::vector<float> &, float &, float &, const std::vector<float> &, const std::vector<float> &);

void downsample(std::vector<float> &, const std::vector<float> &, const int);

void upsample(std::vector<float> &, const std::vector<float> &, const int);

#endif // DY4_FILTER_H
