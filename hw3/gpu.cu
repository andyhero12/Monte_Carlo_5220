#include "common.h"
#include <cuda.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/sort.h>
#include <unistd.h>
#define NUM_THREADS 512

// Put any static global variables here that you will use throughout the simulation.
static int blks;
static constexpr double bin_size =  5.1* cutoff;
// static constexpr double bin_size =  2.1* cutoff;
static int numSMs;
static int grid_len;
__constant__  int grid_len_device;
__constant__  int num_parts_device;
__constant__  double size_device;
static int *d_part_idx;
static int *d_prefix_bins;
static int *d_bin_store_idx; 

__device__ void apply_force_gpu(particle_t& particle, particle_t& neighbor) {
    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;
    if (r2 > cutoff * cutoff)
        return;
    // r2 = fmax( r2, min_r*min_r );
    r2 = (r2 > min_r * min_r) ? r2 : min_r * min_r;
    double r = sqrt(r2);

    //
    //  very simple short-range repulsive force
    //
    double coef = ((1 - cutoff / r) / r2 / mass)*dt;
    atomicAdd(&particle.vx, coef * dx);
    atomicAdd(&particle.vy, coef * dy);
}

__device__ void apply_force_both(particle_t& particle, particle_t& neighbor) {
    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;
    if (r2 > cutoff * cutoff)
        return;
    // r2 = fmax( r2, min_r*min_r );
    r2 = (r2 > min_r * min_r) ? r2 : min_r * min_r;
    double r = sqrt(r2);

    //
    //  very simple short-range repulsive force
    //
    double coef = ((1 - cutoff / r) / r2 / mass)*dt;
    atomicAdd(&particle.vx, coef * dx);
    atomicAdd(&particle.vy, coef * dy);
    atomicAdd(&neighbor.vx, -coef * dx);
    atomicAdd(&neighbor.vy, -coef * dy);
}

__global__ void compute_forces_gpu(particle_t* parts, int* d_part_idx, int* d_prefix_bins) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= num_parts_device)
        return;

        
    particle_t& curP = parts[d_part_idx[tid]];
    int binRow = curP.y / bin_size;
    int binCol = curP.x / bin_size;

    bool goodHigh = binCol != 0;
    bool goodRight = binRow != grid_len_device -1;
    bool goodLow = binCol != grid_len_device -1;
    // curP.ax = curP.ay = 0;
    // double tempvx = 0.0;
    // double tempvy = 0.0;
    int endIndex = d_prefix_bins[(binRow)*grid_len_device + (binCol) +1];
    for (int j = tid+1; j < endIndex; j++){
        apply_force_both(curP, parts[d_part_idx[j]]);
    }
    if (goodLow){
        int startIndex = d_prefix_bins[(binRow)*grid_len_device + (binCol+1)];
        int endIndex = d_prefix_bins[(binRow)*grid_len_device + (binCol+1) +1];
        for (int j = startIndex; j < endIndex; j++){
            apply_force_both(curP, parts[d_part_idx[j]]);
        }
    }
    if (goodRight && goodHigh){
        int startIndex = d_prefix_bins[(binRow+1)*grid_len_device + (binCol-1)];
        int endIndex = d_prefix_bins[(binRow+1)*grid_len_device + (binCol-1) +1];
        for (int j = startIndex; j < endIndex; j++){
            apply_force_both(curP, parts[d_part_idx[j]]);
        }
    }
    if (goodRight){
        int startIndex = d_prefix_bins[(binRow+1)*grid_len_device + (binCol)];
        int endIndex = d_prefix_bins[(binRow+1)*grid_len_device + (binCol) +1];
        for (int j = startIndex; j < endIndex; j++){
            apply_force_both(curP, parts[d_part_idx[j]]);
        }
    }
    if (goodRight && goodLow){
        int startIndex = d_prefix_bins[(binRow+1)*grid_len_device + (binCol+1)];
        int endIndex =   d_prefix_bins[(binRow+1)*grid_len_device + (binCol+1) +1];
        for (int j = startIndex; j < endIndex; j++){
            apply_force_both(curP, parts[d_part_idx[j]]);
        }
    }
}

__global__ void incr_ax(particle_t* particles, int num_parts){
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= 20000000)
        return;
    particle_t* p = &particles[0];
    atomicAdd(&p->ax, 1);
    atomicAdd(&p->ay, 1);
}

__global__ void count_particles_gpu(particle_t* parts, int* d_prefix_bins) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= num_parts_device)
        return;
    particle_t& curP = parts[tid];
    int binRow = curP.y / bin_size;
    int binCol = curP.x / bin_size;
    atomicAdd(d_prefix_bins+(grid_len_device*binRow+binCol),1);
}
__global__ void move_count_particles_gpu(particle_t* particles, int* d_prefix_bins) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= num_parts_device)
        return;
    particle_t* p = &particles[tid];
    p->x += p->vx * dt;
    p->y += p->vy * dt;

    //
    //  bounce from walls
    //
    while (p->x < 0 || p->x > size_device) {
        p->x = p->x < 0 ? -(p->x) : 2 * size_device - p->x;
        p->vx = -(p->vx);
    }
    while (p->y < 0 || p->y > size_device) {
        p->y = p->y < 0 ? -(p->y) : 2 * size_device - p->y;
        p->vy = -(p->vy);
    }
    int binRow = p->y / bin_size;
    int binCol = p->x / bin_size;
    atomicAdd(d_prefix_bins+(grid_len_device*binRow+binCol),1);
    
}

__global__ void place_parts_gpu(particle_t* parts, int* d_bin_store_idx, int* d_part_idx) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= num_parts_device)
        return;
    
    particle_t& curP = parts[tid];
    int binRow = curP.y / bin_size;
    int binCol = curP.x / bin_size;
    int idx = atomicAdd(d_bin_store_idx+grid_len_device*binRow+binCol,1);
    d_part_idx[idx] = tid;
}

inline void first_calc_bins(particle_t* parts){
    cudaMemset(d_prefix_bins, 0, (grid_len*grid_len + 1)*sizeof(int));

    count_particles_gpu<<<blks, NUM_THREADS>>>(parts, d_prefix_bins);

    thrust::exclusive_scan(thrust::device,d_prefix_bins,d_prefix_bins+grid_len*grid_len+1,d_prefix_bins);
    
    cudaMemcpy(d_bin_store_idx, d_prefix_bins, (grid_len*grid_len) * sizeof(int), cudaMemcpyDeviceToDevice);

    place_parts_gpu<<<blks, NUM_THREADS>>>(parts, d_bin_store_idx,d_part_idx);
}

void init_simulation(particle_t* parts, int num_parts, double size) {
    // You can use this space to initialize data objects that you may need
    // This function will be called once before the algorithm begins
    // parts live in GPU memory
    // Do not do any particle simulation here
    cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount,0);
    grid_len = ceil(size/bin_size); 
    cudaMemcpyToSymbol(grid_len_device, &grid_len, sizeof(int));
    cudaMemcpyToSymbol(num_parts_device,&num_parts, sizeof(int));
    cudaMemcpyToSymbol(size_device,&size, sizeof(double));

    cudaMalloc((void**)&d_part_idx,      num_parts * sizeof(int));
    cudaMalloc((void**)&d_prefix_bins,   ((grid_len*grid_len)+1) * sizeof(int));
    cudaMalloc((void**)&d_bin_store_idx, (grid_len*grid_len) * sizeof(int));
    
    blks = (num_parts + NUM_THREADS - 1) / NUM_THREADS;
    incr_ax<<<20000, NUM_THREADS>>>(parts, num_parts);
    first_calc_bins(parts);

}
void simulate_one_step(particle_t* parts, int num_parts, double size) {
    // parts live in GPU memory
    // Rewrite this function
    // Compute forces
    compute_forces_gpu<<<blks, NUM_THREADS>>>(parts, d_part_idx, d_prefix_bins);

    cudaMemset(d_prefix_bins, 0, (grid_len*grid_len + 1)*sizeof(int));

    move_count_particles_gpu<<<blks, NUM_THREADS>>>(parts, d_prefix_bins);

    thrust::exclusive_scan(thrust::device,d_prefix_bins,d_prefix_bins+grid_len*grid_len+1,d_prefix_bins);
    
    cudaMemcpy(d_bin_store_idx, d_prefix_bins, (grid_len*grid_len) * sizeof(int), cudaMemcpyDeviceToDevice);

    place_parts_gpu<<<blks, NUM_THREADS>>>(parts, d_bin_store_idx,d_part_idx);
}
