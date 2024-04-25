#include <chrono>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include "monte-carlo.h"
// Command Line Option Processing
int find_arg_idx(int argc, char** argv, const char* option) {
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], option) == 0) {
            return i;
        }
    }
    return -1;
}

int find_int_arg(int argc, char** argv, const char* option, int default_value) {
    int iplace = find_arg_idx(argc, argv, option);

    if (iplace >= 0 && iplace < argc - 1) {
        return std::stoi(argv[iplace + 1]);
    }

    return default_value;
}
// ==============
// Main Function
// ==============
int main(int argc, char** argv) {
    // Parse Args
    if (find_arg_idx(argc, argv, "-h") >= 0) {
        std::cout << "Options:" << std::endl;
        std::cout << "-s <int>: set num iterations" << std::endl;
        return 0;
    }

    int num_iterations = find_int_arg(argc, argv, "-s", 1000000);
    // Algorithm
    auto start_time = std::chrono::steady_clock::now();
    init_sim();
    result_type result;
#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
    {
        // Calculate the call/put values via Monte Carlo
        monte_carlo_both_price(result,num_iterations);
    }

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

    auto end_time = std::chrono::steady_clock::now();

    std::chrono::duration<double> diff = end_time - start_time;
    double seconds = diff.count();
    // Finalize
    std::cout << "Simulation Time = " << seconds << "\n";
}
