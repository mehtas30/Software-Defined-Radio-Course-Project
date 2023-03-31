/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_FILTER_H
#define DY4_FILTER_H

#include <iostream>
#include <vector>

void impulseResponseLPF(std::vector<float> &, const float, const float, const int, const int);

void impulseResponseBPF(std::vector<float> &, const float, const float, const float, const int);

void resample(std::vector<float> &, std::vector<float> &, const std::vector<float> &, const std::vector<float> &, const int, const int);

void FMDemod(std::vector<float> &, float &, float &, const std::vector<float> &, const std::vector<float> &);

void PLL(std::vector<float> &, const float, const float, const float, const float, const float, float &, float &, float &, float &, float &, float &);

void mixer(std::vector<float> &, const std::vector<float> &, const std::vector<float> &);

void LRExtraction(std::vector<float> &, std::vector<float> &, const std::vector<float> &, const std::vector<float> &);

#endif // DY4_FILTER_H
