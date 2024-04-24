//=========================================================================
// monte-carlo_test.cpp
//=========================================================================
// @brief: testbench for Monte Carlo Xcel

#include <iostream>
#include <cmath>

#include "monte-carlo.h"

int main(int argc, char **argv)
{
  // Calculate the call/put values via Monte Carlo
  result_type result;
  monte_carlo_both_price(result);
  theta_type call = result.call;
  theta_type put = result.put;
  // Finally we output the parameters and prices
  std::cout << "Testbench" << std::endl;
  std::cout << "Number of Paths: " << num_sims << std::endl;
  std::cout << "Underlying:      " << S << std::endl;
  std::cout << "Strike:          " << K << std::endl;
  std::cout << "Risk-Free Rate:  " << r << std::endl;
  std::cout << "Volatility:      " << v << std::endl;
  std::cout << "Maturity:        " << T << std::endl;

  std::cout << "Call Price:      " << call << std::endl;
  std::cout << "Put Price:       " << put << std::endl;

  return 0;
}
