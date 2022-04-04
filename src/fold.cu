#include <exception>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <random>
#include <cstdio>
#include <iostream>
#include <Eigen/Dense>
#include <curand_kernel.h>
#include <curand.h>
#include <mpi.h>
#include <cuda_runtime_api.h>
#include <cuda.h>


#include "protein.hpp"

__device__ bool check_residue(Protein protein, Eigen::Vector3i pos) {
    for (int i = 0; i < protein.sz; i++) {
        if (protein.positions[i] == pos) return false;
    }
    return true;
}
__device__ void step(Protein& protein, curandState* state) {
    auto index = curand(state) % protein.sz;
    auto pos = protein.positions[index];
    Eigen::Vector3i moves[8];
    size_t len_moves = 0;
    if (index == 0 || index == protein.sz) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 2; i++) {
                Eigen::Vector3i new_pos = pos;
                (j == 0) ? new_pos[i] += 1 : new_pos[i] -= 1;
                if (check_residue(protein, new_pos) == true) {
                    moves[len_moves+i] = new_pos;
                    len_moves++;
                }
            }
        }
    }
}

__device__ int get_cost(Protein protein) {
    int cost = 0;
    const Eigen::Vector3i unit[3] = {{0,0,1}, {0,1,0}, {1,0,0}};
    for (int i = 0; i < protein.sz; i++) {
        for (int i = 0; i < 3; i++) {
            if (check_residue(protein, protein.positions[i] + unit[i])) cost++;
            if (check_residue(protein, protein.positions[i] - unit[i])) cost++;
        }
    }
    return cost;
}

__global__ void __anneal_multistart_singlestrat(Sequence seq, Protein* proteins, int* min_cost, curandState* states) {
    // assert(get_cost(*(proteins+threadIdx.x)));
    int cost = get_cost(*(proteins+threadIdx.x));
    // assert(cost == 239);
    // assert(0);
    for (int i = 0; i < 1000; i++) {
        step(proteins[threadIdx.x], &states[threadIdx.x]);
    }
    cost = get_cost(*(proteins+threadIdx.x));
    for (int offset = 32; offset > 0; offset /= 2) {
        auto other = __shfl_down_sync(0xFFFFFFFF, cost, offset);
        cost = cost > other ? cost : other;
    }
    // assert(cost == 3289);
    if (threadIdx.x == 0) *min_cost = cost;
}

// void copy_protein(Protein& p, void* buffer) {
//      cudaMallocManaged(&p.positions, p.sequence.size*sizeof(Eigen::Vector3i));
//      cudaMallocManaged(&p.sequence, p.sequence.size*sizeof(Amino));
// }

std::tuple<std::vector<Protein>, int> anneal_multistart_singlestrat(Sequence seq) {
    Protein* dev_proteins;
    cudaMallocManaged(&dev_proteins, sizeof(Protein) * 32);
    for (int i = 0; i < 32; i++) dev_proteins[i] = Protein::random(seq);
    curandState *dev_random;
    cudaMalloc(&dev_random, 32*sizeof(curandState));	
    int* cost;
    cudaMalloc(&cost, sizeof(int));
    __anneal_multistart_singlestrat<<<1,32>>>(seq, dev_proteins, cost, dev_random);
    cudaDeviceSynchronize();
    int local_cost;
    cudaMemcpy(&local_cost, cost, 4, cudaMemcpyDeviceToHost);
    std::vector<Protein> out;
    for (int i = 0; i < 32; i++) {
	out.push_back(dev_proteins[i]);
    }
    return std::make_tuple(out, local_cost);
}
