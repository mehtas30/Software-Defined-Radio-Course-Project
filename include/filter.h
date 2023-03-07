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
void impulseResponseLPF(float, float, unsigned short int, std::vector<float> &);
void convolveFIR(std::vector<float> &, const std::vector<float> &, const std::vector<float> &);
void demodFM(const std::vector<float> &i_ds, const std::vector<float> &q_ds, std::vector<float> &demod, float p_i=0, float p_q=0);

#endif // DY4_FILTER_H
