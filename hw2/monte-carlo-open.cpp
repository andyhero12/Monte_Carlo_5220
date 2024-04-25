#include <iostream>
#include "monte-carlo.h"
#include <random>
#include <omp.h>
std::default_random_engine generator;
std::uniform_real_distribution<double> distribution(0.0,1.0);
// Function to generate thread-local random number generators
std::default_random_engine create_local_generator()
{
    std::random_device rd;
    std::default_random_engine generator(rd());
    return generator;
}
// Box muller algorithm
theta_type gaussian_box_muller(std::default_random_engine& generator)
{
    theta_type x = 0.45543;
    theta_type euclid_sq = 0.353308;
    theta_type euclid_sq_temp = 0.0;
    theta_type temp_x = 0;
    theta_type temp_y = 0;

    // Continue generating two uniform random variables
    // until the square of their "euclidean distance"
    // is less than unity
    for (int i = 0; i < 20; i++)
    {
        temp_x = distribution(generator) * 2.0 - 1;
        temp_y = distribution(generator) * 2.0 - 1;
        euclid_sq_temp = temp_x * temp_x + temp_y * temp_y;
        if (euclid_sq_temp < 1.0)
        {
            euclid_sq = euclid_sq_temp;
            x = temp_x;
        }
    }

    return x * std::sqrt(-2 * std::log(euclid_sq) / euclid_sq);
}

void init_sim(){

}

// Pricing a European vanilla option with a Monte Carlo method
void monte_carlo_both_price(result_type &result, long long iterations){
  const theta_type S_adjust = S * std::exp(T * (r - 0.5 * v * v));
  const theta_type sqrt_const = std::sqrt(v * v * T);
  theta_type call_payoff_sum = 0.0;
  theta_type put_payoff_sum = 0.0;
  omp_set_num_threads(64);
#pragma omp parallel reduction(+: call_payoff_sum, put_payoff_sum)
    {
        std::default_random_engine generator = create_local_generator();
#pragma omp for
        for (long long i = 0; i < iterations; i++)
        {
            theta_type gauss_bm = gaussian_box_muller(generator);
            theta_type S_cur = S_adjust * std::exp(sqrt_const * gauss_bm);
            theta_type zero1 = 0.0;
            theta_type zero2 = 0.0;
            theta_type call_val = S_cur - K;
            theta_type put_val = K - S_cur;
            call_payoff_sum += std::max(call_val, zero1);
            put_payoff_sum += std::max(put_val, zero2);
        }
    }

    theta_type cast_num_sims = iterations;
    theta_type call = (call_payoff_sum / cast_num_sims) * std::exp(-r * T);
    theta_type put = (put_payoff_sum / cast_num_sims) * std::exp(-r * T);

    result.call = call;
    result.put = put;
}