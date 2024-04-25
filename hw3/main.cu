#include <chrono>
#include <cmath>
#include <cstring>
#include <cuda.h>
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


    result_type* callPut = new result_type;
    callPut->call = 0.0;
    callPut->put = 0.0;

    result_type* callPut_gpu;
    cudaMalloc((void**)&callPut_gpu, sizeof(result_type));
    cudaMemcpy(callPut_gpu, callPut, sizeof(result_type), cudaMemcpyHostToDevice);

    // Algorithm
    auto start_time = std::chrono::steady_clock::now();

    init_sim();

    // cudaMemcpy(callPut, callPut_gpu, sizeof(result_type), cudaMemcpyDeviceToHost);
    monte_carlo_both_price(callPut, num_iterations);
    cudaDeviceSynchronize();
    auto end_time = std::chrono::steady_clock::now();

    std::chrono::duration<double> diff = end_time - start_time;
    double seconds = diff.count();

    theta_type call = callPut->call;
    theta_type put = callPut->put;
    // Finally we output the parameters and prices
    std::cout << "Testbench" << std::endl;
    std::cout << "Number of Paths: " << num_iterations  << std::endl;
    std::cout << "Underlying:      " << S << std::endl;
    std::cout << "Strike:          " << K << std::endl;
    std::cout << "Risk-Free Rate:  " << r << std::endl;
    std::cout << "Volatility:      " << v << std::endl;
    std::cout << "Maturity:        " << T << std::endl;

    std::cout << "Call Price:      " << call << std::endl;
    std::cout << "Put Price:       " << put << std::endl;
    // Finalize
    std::cout << "Simulation Time = " << seconds << " seconds for " << num_iterations << " iterations.\n";
    cudaFree(callPut_gpu);
    delete[] callPut;
}
