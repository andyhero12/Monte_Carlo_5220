#include <iostream>
#include "monte-carlo.h"

constexpr theta_type rand_two_div_max = 2.0 / RAND_MAX;
// Function to generate a random number in the range [0, 1)
theta_type generate_rand() {
    // theta_type casted_seed = pseudo_random(lfsr1);
    theta_type casted_seed = rand();
    return rand_two_div_max * casted_seed - 1;
}


// Box muller algorithm
theta_type gaussian_box_muller()
{
  theta_type x = 0.45543;
  theta_type y = -0.337388;
  theta_type euclid_sq = 0.353308;
  
  theta_type euclid_sq_temp = 0.0;
  theta_type epsilon = 0.00001;
  theta_type temp_x = 0;
  theta_type temp_y = 0;

  // Continue generating two uniform random variables
  // until the square of their "euclidean distance"
  // is less than unity
  for (int i = 0; i < 20; i++) {
    temp_x = generate_rand();
    temp_y = generate_rand();
    euclid_sq_temp = temp_x * temp_x + temp_y * temp_y;
    if (euclid_sq_temp < 1.0) {
      euclid_sq = euclid_sq_temp;
      x = temp_x;
      y = temp_y;
    }
  }
  
  return x * sqrt(-2 * log(euclid_sq) / euclid_sq);
}

void init_sim(){
}
// Pricing a European vanilla option with a Monte Carlo method
void monte_carlo_both_price(result_type &result, int iterations)
{

  const theta_type S_adjust = S * exp(T * (r - 0.5 * v * v));
  const theta_type sqrt_const = sqrt(v * v * T);
  theta_type S_cur = 0.0;
  theta_type call_payoff_sum = 0.0;
  theta_type put_payoff_sum = 0.0;
  
  for (int i = 0; i < iterations; i++) {
    theta_type gauss_bm = gaussian_box_muller();
    S_cur = S_adjust * exp(sqrt_const * gauss_bm);
    theta_type zero1 = 0.0;
    theta_type zero2 = 0.0;
    theta_type call_val = S_cur - K;
    theta_type put_val = K - S_cur;
    call_payoff_sum += fmax(call_val, zero1);
    put_payoff_sum += fmax(put_val, zero2);
  }
  

  theta_type cast_num_sims = num_sims;
  theta_type call = (call_payoff_sum / cast_num_sims) * exp(-r * T);
  theta_type put = (put_payoff_sum / cast_num_sims) * exp(-r * T);

  result.call = call;
  result.put = put;
}
