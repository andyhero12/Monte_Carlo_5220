#ifndef MONTE_H
#define MONTE_H

#include <cmath>
#include <iostream>

typedef struct result_type result_type;

typedef float theta_type;

static constexpr int num_sims = 1000000;      // Number of simulated asset paths
static constexpr theta_type S = 100.0;        // Option price
static constexpr theta_type K = 100.0;        // Strike price
static constexpr theta_type r = 0.05;         // Risk-free rate (5%)
static constexpr theta_type v = 0.2;          // Volatility of the underlying (20%)
static constexpr theta_type T = 1.0;          // One year until expiry


theta_type gaussian_box_muller();

void init_sim();
void monte_carlo_both_price(result_type &result, int iterations);

struct result_type {
theta_type call; 
theta_type put;
};


#endif // MONTE_H
