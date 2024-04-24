#include <iostream>
#include "hls_math.h"
#include "monte-carlo.h"

int num_sims = 1000000;      // Number of simulated asset paths
theta_type S = 100.0;        // Option price
theta_type K = 100.0;        // Strike price
theta_type r = 0.05;         // Risk-free rate (5%)
theta_type v = 0.2;          // Volatility of the underlying (20%)
theta_type T = 1.0;          // One year until expiry

void dut(hls::stream<bit32_t> &strm_in, hls::stream<bit32_t> &strm_out) {
  // Read the results from the module
  int dummy = strm_in.read();

  result_type result;
  monte_carlo_both_price(result);
  theta_type call = result.call;
  theta_type put = result.put;

  union { float fval; int ival; } ucall;
  union { float fval; int ival; } uput;

  ucall.fval = call;
  uput.fval = put;

  int icall = ucall.ival;
  int iput = uput.ival;

  // write output to stream
  strm_out.write(icall);
  strm_out.write(iput);
}

ap_uint<32> lfsr1 = 0xdcba;  // Initial seed
ap_uint<32> lfsr2 = 1234;  // Initial seed

// PRNG - pseduo random number generator using linear feedback shift registers
unsigned int pseudo_random(ap_uint<32>& lfsr) {
  bool b_32 = lfsr.get_bit(32-32);
  bool b_22 = lfsr.get_bit(32-22);
  bool b_2 = lfsr.get_bit(32-2);
  bool b_1 = lfsr.get_bit(32-1);
  bool new_bit = b_32 ^ b_22 ^ b_2 ^ b_1;
  lfsr = lfsr >> 1;
  lfsr.set_bit(31, new_bit);
  return lfsr.to_uint();
}

constexpr theta_type rand_two_div_max = 2.0 / RAND_MAX;
// Function to generate a random number in the range [0, 1)
theta_type generate_rand1() {
    theta_type casted_seed = pseudo_random(lfsr1);
    return rand_two_div_max * casted_seed - 1;
}

// Function to generate a random number in the range [0, 1)
theta_type generate_rand2() {
    theta_type casted_seed = pseudo_random(lfsr2);
    return rand_two_div_max * casted_seed - 1;
}

// custom log function 
template <typename T>
T custom_log(const T& x)
{
  if (x <= 0)
  {
    std::cerr << "Error: Input must be greater than 0" << std::endl;
    return -1.0; // Error value
  }
  
  const int logTerms = 10;

  T result = 0.0;
  T term = (x - 1) / (x + 1);
  T term_squared = term * term;
  T numerator = term;
  T denominator = 1;

  LOG_LABEL:
  for (int i = 1; i <= logTerms; i++)
  {
    result += numerator / denominator;
    numerator *= term_squared;
    denominator += 2;
  }

  return 2 * result;
}

// custom exp function
template <typename T>
T custom_exp(const T& x)
{
  T result = 1.0;
  T term = 1.0;
  T factorial = 1.0;
  const int expTerms = 10;

  EXP_LABEL:
  for (int i = 1; i <= expTerms; i++)
  {
    term *= x / i;
    result += term;
  }

  return result;
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

  GAUSS_LABEL:
  for (int i = 0; i < 20; i++) {
    temp_x = generate_rand1();
    temp_y = generate_rand2();
    euclid_sq_temp = temp_x * temp_x + temp_y * temp_y;
    if (euclid_sq_temp < 1.0) {
      euclid_sq = euclid_sq_temp;
      x = temp_x;
      y = temp_y;
    }
  }
  
  return x * hls::sqrt(-2 * custom_log<theta_type>(euclid_sq) / euclid_sq);
}

// Pricing a European vanilla option with a Monte Carlo method
void monte_carlo_both_price(result_type &result)
{

  const theta_type S_adjust = S * custom_exp<theta_type>(T * (r - 0.5 * v * v));
  const theta_type sqrt_const = hls::sqrt(v * v * T);
  theta_type S_cur = 0.0;
  theta_type call_payoff_sum = 0.0;
  theta_type put_payoff_sum = 0.0;
  
  GAUSS_GEN_LABEL:
  for (int i = 0; i < 1000000; i++) {
    theta_type gauss_bm = gaussian_box_muller();
    S_cur = S_adjust * custom_exp<theta_type>(sqrt_const * gauss_bm);
    theta_type zero1 = 0.0;
    theta_type zero2 = 0.0;
    theta_type call_val = S_cur - K;
    theta_type put_val = K - S_cur;
    call_payoff_sum += hls::fmax(call_val, zero1);
    put_payoff_sum += hls::fmax(put_val, zero2);
  }
  

  theta_type cast_num_sims = num_sims;
  theta_type call = (call_payoff_sum / cast_num_sims) * custom_exp<theta_type>(-r * T);
  theta_type put = (put_payoff_sum / cast_num_sims) * custom_exp<theta_type>(-r * T);

  result.call = call;
  result.put = put;
}
