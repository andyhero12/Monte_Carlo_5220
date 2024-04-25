#include "common.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#pragma GCC optimize ("O3")
#pragma GCC target ("arch=znver3")

typedef std::vector<particle_t*> bin;

size_t grid_len;
std::vector<bin> bins; 
void org(particle_t* parts, int num_parts, double size);

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
    particle.ax += coef * dx;
    particle.ay += coef * dy;

    neighbor.ax -= coef * dx;
    neighbor.ay -= coef * dy;

}

static inline void clear_bins(std::vector<bin>& curBins){
    for (bin& cur_bin : curBins){
        cur_bin.clear();
    }
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


// grid_len x grid_len
// row major
// 0,  1,  2,  3, 
// 4,  5,  6,  7, 
// 8,  9,  10, 11, 
// 12, 13, 14, 15
void reset_accel(particle_t* parts, int num_parts){
    for (int i =0 ; i< num_parts;i++){
        particle_t& curP = parts[i];
        curP.ax = curP.ay = 0;
    }
}
void reset_put_in_bins(particle_t* parts, int num_parts){
    // do logic to put into bin
    // q = (x + y - 1) / y;
    clear_bins(bins);
    for (int i = 0; i< num_parts; i++){
        particle_t& curP = parts[i];
        // reset
        int binRow = curP.y / cutoff;
        int binCol = curP.x / cutoff;
        bins[binRow * grid_len + binCol].push_back(parts+i);
    }
}
void move_reset_put_in_bins(particle_t* parts, int num_parts, double size){
    // do logic to put into bin
    // q = (x + y - 1) / y;
    if (num_parts <= 100000){
        clear_bins(bins);
        for (int i = 0; i< num_parts; i++){
            particle_t& curP = parts[i];
            // move
            move(curP, size);
            // reset
            curP.ax = curP.ay = 0;

            int binRow = curP.y / cutoff;
            int binCol = curP.x / cutoff;
            bins[binRow * grid_len + binCol].emplace_back(parts + i);
        } 
    }else{
        for (int i = 0; i< num_parts; i++){
            particle_t& curP = parts[i];
            int binRowOrg = curP.y / cutoff;
            int binColOrg = curP.x / cutoff;
            int orgInd = binRowOrg*grid_len + binColOrg;
            move(curP, size);
            curP.ax = curP.ay = 0;
            int binRow = curP.y / cutoff;
            int binCol = curP.x / cutoff;
            int ind = binRow * grid_len + binCol;
            if (orgInd != ind){
                auto pos = std::find(bins[orgInd].begin(), bins[orgInd].end(), &curP);
                bins[orgInd].erase(pos);
                bins[ind].emplace_back(&curP);
            }
        }
    }
}

static inline void compare_bins(const bin& orgBin, const bin& neighBin){
    for (particle_t* org: orgBin){
        for (particle_t* dest: neighBin){
            apply_force_both(*org, *dest);
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
static inline void bin_simulate_one_step(particle_t* parts, int num_parts, double size) {
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
    move_reset_put_in_bins(parts, num_parts,size);
}

void init_simulation(particle_t* parts, int num_parts, double size) {
    // You can use this space to initialize static, global data objects
    // that you may need. This function will be called once before the
    // algorithm begins. Do not do any particle simulation here
    grid_len = ceil(size/cutoff);
    bins.resize(grid_len*grid_len);
    reset_accel(parts, num_parts);  
    reset_put_in_bins(parts, num_parts);
}

void simulate_one_step(particle_t* parts, int num_parts, double size) {
    bin_simulate_one_step(parts, num_parts,size);
    // org(parts, num_parts, size);
}

void org(particle_t* parts, int num_parts, double size){
    for (int i = 0; i < num_parts; ++i) {
        for (int j = i+1; j < num_parts; ++j) {
            apply_force_both(parts[i], parts[j]);
        }
    }

    // Move Particles
    for (int i = 0; i < num_parts; ++i) {
        move(parts[i], size);
        parts[i].ax = parts[i].ay = 0;
    }
}