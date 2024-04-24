#pragma once

#ifndef MONTE_H
#define MONTE_H

#include <cmath>
#include <iostream>

#include "typedefs.h"
#include <hls_stream.h>

typedef struct result_type result_type;

typedef float theta_type;
extern int num_sims;
extern theta_type S, K, r, v, T;


void dut(hls::stream<bit32_t> &strm_out);

theta_type generate_rand();

theta_type gaussian_box_muller();

void monte_carlo_both_price(result_type &result);

struct result_type {
theta_type call; 
theta_type put;
};


#endif // MONTE_H
