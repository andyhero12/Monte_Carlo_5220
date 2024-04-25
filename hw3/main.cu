#include "common.h"
#include <chrono>
#include <cmath>
#include <cstring>
#include <cuda.h>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include "monte-carlo.h"

// Particle Initialization
void init_particles(particle_t* parts, int num_parts, double size, int part_seed) {
    std::random_device rd;
    std::mt19937 gen(part_seed ? part_seed : rd());

    int sx = (int)ceil(sqrt((double)num_parts));
    int sy = (num_parts + sx - 1) / sx;

    std::vector<int> shuffle(num_parts);
    for (int i = 0; i < shuffle.size(); ++i) {
        shuffle[i] = i;
    }

    for (int i = 0; i < num_parts; ++i) {
        // Make sure particles are not spatially sorted
        std::uniform_int_distribution<int> rand_int(0, num_parts - i - 1);
        int j = rand_int(gen);
        int k = shuffle[j];
        shuffle[j] = shuffle[num_parts - i - 1];

        // Distribute particles evenly to ensure proper spacing
        parts[i].x = size * (1. + (k % sx)) / (1 + sx);
        parts[i].y = size * (1. + (k / sx)) / (1 + sy);

        // Assign random velocities within a bound
        std::uniform_real_distribution<float> rand_real(-1.0, 1.0);
        parts[i].vx = rand_real(gen);
        parts[i].vy = rand_real(gen);
    }
}

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
