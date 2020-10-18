#ifndef RIR_GENERATOR_CORE_H
#define RIR_GENERATOR_CORE_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

void computeRIR(double* imp, double c, double fs, double* rr, int nr_of_mics, int nsamples, double* ss, double* LL, double* beta, char mtype, int order, double* angle, int hp_filter);

#ifdef __cplusplus
}
#endif

#endif
