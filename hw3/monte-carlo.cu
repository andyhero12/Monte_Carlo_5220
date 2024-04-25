#include <iostream>
#include "monte-carlo.h"
#include <curand_kernel.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/sort.h>
#include <random>
#define NUM_THREADS 256
static int numSMs;
static theta_type* calls;
static theta_type* putOptionValues;
std::default_random_engine generator;
std::uniform_real_distribution<double> distribution(0.0,1.0);

void init_sim(){
  cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount,0);
  cudaMalloc((void**)&calls, 32 * numSMs *  sizeof(theta_type));
  cudaMemset(calls, 0, 32 * numSMs* sizeof(theta_type));
  cudaMalloc((void**)&putOptionValues, 32 * numSMs *  sizeof(theta_type));
  cudaMemset(putOptionValues, 0, 32 * numSMs* sizeof(theta_type));

}

// Box muller algorithm
__device__
static inline theta_type device_box_muller(curandState* state)
{
  theta_type x = 0.45543;
  theta_type euclid_sq = 0.353308;
  theta_type euclid_sq_temp = 0.0;
  theta_type temp_x = 0;
  theta_type temp_y = 0;

  // Continue generating two uniform random variables
  // until the square of their "euclidean distance"
  // is less than unity
  for (int i = 0; i < 20; i++) {
    temp_x = curand_uniform(state) * 2.0 - 1;
    temp_y = curand_uniform(state) * 2.0 - 1;
    euclid_sq_temp = temp_x * temp_x + temp_y * temp_y;
    if (euclid_sq_temp < 1.0) {
      euclid_sq = euclid_sq_temp;
      x = temp_x;
    }
  }
  
  return x * sqrt(-2 * log(euclid_sq) / euclid_sq);
}
__global__ 
void monte_carlo_iteration(theta_type* calls, theta_type* putOption, int iterations, int seed){

  __shared__ theta_type callValues[NUM_THREADS];
  __shared__ theta_type putValues[NUM_THREADS];
  if (threadIdx.x == 0){
    for (int i =0 ; i< NUM_THREADS; i++){
      putValues[i] = 0.0;
      callValues[i] = 0.0;
    }
  }
  __syncthreads();
  int index = blockDim.x * blockIdx.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  const theta_type S_adjust = S * exp(T * (r - 0.5 * v * v));
  const theta_type sqrt_const = sqrt(v * v * T);
  theta_type call_payoff_sum = 0.0;
  theta_type put_payoff_sum = 0.0;
  curandState state;
  curand_init(seed + blockIdx.x, 0, 0, &state);
  for (int i = index; i< iterations; i += stride){
    theta_type gauss_bm = device_box_muller(&state);
    theta_type S_cur = S_adjust * exp(sqrt_const * gauss_bm);
    theta_type zero1 = 0.0;
    theta_type zero2 = 0.0;
    theta_type call_val = S_cur - K;
    theta_type put_val = K - S_cur;
    call_payoff_sum += max(call_val, zero1);
    put_payoff_sum += max(put_val, zero2);
    // Perform the Monte Carlo simulation for this iteration
  }
  callValues[threadIdx.x] = call_payoff_sum;
  putValues[threadIdx.x] = put_payoff_sum;
  __syncthreads();
  unsigned int tid = threadIdx.x;
  theta_type builderCall = 0.0;
  theta_type builderPut = 0.0;

  if (tid == 0){
    for (int i =0; i< NUM_THREADS; i++){
      builderCall += callValues[i];
      builderPut += putValues[i];
    }
    calls[blockIdx.x] = builderCall;
    putOption[blockIdx.x] = builderPut;
  }
}

void monte_carlo_both_price(result_type* result, int iterations)
{
  monte_carlo_iteration<<<32 * numSMs, NUM_THREADS>>>(calls, putOptionValues, iterations, time(NULL));
  
  theta_type total_call_sum = 0.0;
  theta_type total_put_sum = 0.0;
  total_call_sum = thrust::reduce(thrust::device, calls, calls + 32 * numSMs,0.0);
  total_put_sum = thrust::reduce(thrust::device, putOptionValues, putOptionValues + 32 * numSMs,0.0);
  theta_type cast_num_sims = iterations;
  theta_type call = (total_call_sum / cast_num_sims) * exp(-r * T);
  theta_type put = (total_put_sum / cast_num_sims) * exp(-r * T);

  result->call = call;
  result->put = put;
}
