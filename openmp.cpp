#include "common.h"
#include <omp.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <utility>
#pragma GCC optimize ("O3")
#pragma GCC target ("arch=znver3")

#define NUM_THREADS 64
// Put any static global variables here that you will use throughout the simulation.
// typedef std::vector<particle_t*> bin ;
typedef std::vector<particle_t*> bin;
size_t grid_len __attribute__((aligned(64)));

// bin* bins __attribute__((aligned(64))); 
// bin* nxtBins __attribute__((aligned(64))) ;
// omp_lock_t* locks __attribute__((aligned(64)));
std::vector<std::pair<double,double>> accelTemp;
std::vector<bin> bins __attribute__((aligned(64))); 
std::vector<bin> nxtBins __attribute__((aligned(64))) ;
std::vector<omp_lock_t> locks __attribute__((aligned(64)));


// Apply the force from neighbor to particle
static inline void apply_force(particle_t& particle, particle_t& neighbor) {
    // Calculate Distance
    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;

    // Check if the two particles should interact
    if (r2 > cutoff * cutoff)
        return;

    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);

    // Very simple short-range repulsive force
    double coef = (1 - cutoff / r) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
}

static inline void apply_force_both(particle_t& particle, particle_t& neighbor) {
    // Calculate Distance
    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;

    // Check if the two particles should interact
    if (r2 > cutoff * cutoff)
        return;

    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);

    // Very simple short-range repulsive force
    double coef = (1 - cutoff / r) / r2 / mass;
    #pragma omp atomic
    particle.ax += coef * dx;
    #pragma omp atomic
    particle.ay += coef * dy;
    #pragma omp atomic
    neighbor.ax -= coef * dx;
    #pragma omp atomic
    neighbor.ay -= coef * dy;

}
// Integrate the ODE
static inline void move(particle_t& p, double size) {
    // Slightly simplified Velocity Verlet integration
    // Conserves energy better than explicit Euler method
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x += p.vx * dt;
    p.y += p.vy * dt;

    // Bounce from walls
    while (p.x < 0 || p.x > size) {
        p.x = p.x < 0 ? -p.x : 2 * size - p.x;
        p.vx = -p.vx;
    }

    while (p.y < 0 || p.y > size) {
        p.y = p.y < 0 ? -p.y : 2 * size - p.y;
        p.vy = -p.vy;
    }
}

void reset_put_in_bins(particle_t* parts, int num_parts){
    for (int i = 0; i< num_parts; i++){
        particle_t& curP = parts[i];
        // reset
        curP.ax = curP.ay = 0;

        int binRow = curP.y / cutoff;
        int binCol = curP.x / cutoff;
        bins[binRow * grid_len + binCol].emplace_back(&curP);
    }
}

static inline void move_reset_put_in_bins(particle_t* parts, int num_parts, double size){
    #pragma omp for
    for (int i = 0; i< num_parts; i++){
        particle_t& curP = parts[i];
        int binRowOrg = curP.y / cutoff;
        int binColOrg = curP.x / cutoff;
        int orgInd = binRowOrg*grid_len + binColOrg;
        move(curP, size);
        // reset
        curP.ax = curP.ay = 0;
        int binRow = curP.y / cutoff;
        int binCol = curP.x / cutoff;
        int ind = binRow * grid_len + binCol;
        if (orgInd != ind){
            omp_set_lock(&locks[orgInd]);
            auto pos = std::find(bins[orgInd].begin(), bins[orgInd].end(), &curP);
            bins[orgInd].erase(pos);
            omp_unset_lock(&locks[orgInd]);
            omp_set_lock(&locks[ind]);
            bins[ind].emplace_back(&curP);
            omp_unset_lock(&locks[ind]);
        }
    }
    #pragma omp barrier
}

static inline void compare_bins(const bin& orgBin, const bin& neighBin){
    for (particle_t* org: orgBin){
        for (particle_t* dest: neighBin){
            apply_force_both(*org, *dest);
        }
    }
}
static inline void compare_bins_one(const bin& orgBin, const bin& neighBin){
    for (particle_t* org: orgBin){
        for (particle_t* dest: neighBin){
            apply_force(*org, *dest);
        }
    }
}
static inline void compare_bins_same(const bin& orgBin){
    for (int i = 0; i< orgBin.size();i++){
        for (int j = i + 1; j< orgBin.size();j++){
            apply_force_both(*orgBin[i],*orgBin[j]);
        }
    }
}

void init_simulation(particle_t* parts, int num_parts, double size) {
	// You can use this space to initialize data objects that you may need
	// This function will be called once before the algorithm begins
	// Do not do any particle simulation here
    // omp_set_num_threads(NUM_THREADS);
    grid_len = ceil(size/cutoff);
    bins.resize(grid_len*grid_len);
    nxtBins.resize(grid_len*grid_len);
    locks.resize(grid_len*grid_len);
    accelTemp.resize(num_parts);
    // bins    = (bin*) aligned_alloc(64, sizeof(bin) * grid_len * grid_len);
    // nxtBins = (bin*) aligned_alloc(64, sizeof(bin) * grid_len * grid_len);
    // locks   = (omp_lock_t*) aligned_alloc(64, sizeof(omp_lock_t) * grid_len * grid_len);
    for (int i =0 ; i < grid_len * grid_len; i++){
        // bins[i] = bin();
        // nxtBins[i] = bin();
        // locks[i] = omp_lock_t();
        omp_init_lock(&locks[i]);
    }
    reset_put_in_bins(parts, num_parts);
}

static inline void bin_simulate_one_step(particle_t* parts, int num_parts, double size) {
    #pragma omp for
    for (int r = 0 ; r < grid_len; r++){
        for (int c = 0 ;c < grid_len; c++){
            const bin& cur_bin = bins[c + r*grid_len];
            // top left bin
            if (r != 0 && c != 0 ){
                const bin& neighBin = bins[(c-1) + (r-1)*grid_len];
                compare_bins_one(cur_bin,neighBin);
            }
            // top bin
            if (r != 0 ){
                const bin& neighBin = bins[c + (r-1)*grid_len];
                compare_bins_one(cur_bin,neighBin);
            }
            // top right bin
            if (r != 0 && c != grid_len - 1 ){
                const bin& neighBin = bins[(c+1) + (r-1)*grid_len];
                compare_bins_one(cur_bin,neighBin);
            }
            // left bin
            if (c != 0){
                const bin& neighBin = bins[(c-1) + r*grid_len];
                compare_bins_one(cur_bin,neighBin);
            }
            // right bin
            if (c != grid_len - 1){
                const bin& neighBin = bins[(c+1) + r*grid_len];
                compare_bins_one(cur_bin,neighBin);
            }
            // below left bin
            if (r != grid_len - 1 && c != 0){
                const bin& neighBin = bins[(c-1) + (r+1)*grid_len];
                compare_bins_one(cur_bin,neighBin);
            }
            // below bin
            if (r != grid_len - 1){
                const bin& neighBin = bins[c + (r+1)*grid_len];
                compare_bins_one(cur_bin,neighBin);
            }
            // below right bin
            if (r != grid_len - 1 && c != grid_len - 1){
                const bin& neighBin = bins[(c+1) + (r+1)*grid_len];
                compare_bins_one(cur_bin,neighBin);
            }
            // same bin
            compare_bins_one(cur_bin,cur_bin);
        }
    }
}


static inline void working(particle_t* parts, int num_parts, double size){
    // auto org_start_time  = std::chrono::steady_clock::now();
    #pragma omp for
    for (int r = 0 ; r < grid_len; r++){
        for (int c = 0 ;c < grid_len; c++){
            const bin& cur_bin = bins[c + r*grid_len];
            // right bin
            if (c != grid_len - 1){
                const bin& neighBin = bins[(c+1) + r*grid_len];
                compare_bins(cur_bin,neighBin);
            }
            // below left bin
            if (r != grid_len - 1 && c != 0){
                const bin& neighBin = bins[(c-1) + (r+1)*grid_len];
                compare_bins(cur_bin,neighBin);
            }
            // below bin
            if (r != grid_len - 1){
                const bin& neighBin = bins[c + (r+1)*grid_len];
                compare_bins(cur_bin,neighBin);
            }
            // below right bin
            if (r != grid_len - 1 && c != grid_len - 1){
                const bin& neighBin = bins[(c+1) + (r+1)*grid_len];
                compare_bins(cur_bin,neighBin);
            }
            // same bin
            compare_bins_same(cur_bin);
        }
    }
    #pragma omp barrier
    // auto start_time  = std::chrono::steady_clock::now();
    // std::chrono::duration<double> diff2 = start_time - org_start_time;
    // double seconds2 = diff2.count();
    // #pragma omp master
    // {
    //     printf("computeTime : %f\n", seconds2);
    // }

    move_reset_put_in_bins(parts, num_parts,size);

    // auto end_time = std::chrono::steady_clock::now();
    // std::chrono::duration<double> diff = end_time - start_time;
    // double seconds = diff.count();

    // #pragma omp master
    // {
    //     printf("movetime : %f\n", seconds);
    // }
}
void simulate_one_step(particle_t* parts, int num_parts, double size) {
    // org(parts,num_parts, size);  
    working(parts, num_parts, size);
    // bin_simulate_one_step(parts, num_parts, size);
}

