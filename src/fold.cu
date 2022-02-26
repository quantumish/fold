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

#include "protein.hpp"

__device__ void step(Protein& protein, curandState* state) {    
    // rand() % 
}

__device__ bool check_residue(Protein protein, Eigen::Vector3i pos) {
    for (int i = 0; i < protein.sequence.size; i++) {
		if (protein.positions[i] == pos) return false;
    }
    return true;
}

__device__ int get_cost(Protein protein) {
    int cost = 0;
    const Eigen::Vector3i unit[3] = {{0,0,1}, {0,1,0}, {1,0,0}};
    for (int i = 0; i < protein.sequence.size; i++) {
		for (int i = 0; i < 3; i++) {
			if (check_residue(protein, protein.positions[i] + unit[i])) cost++;
			if (check_residue(protein, protein.positions[i] - unit[i])) cost++;	    
		}
    }
}

__global__ void __anneal_multistart_singlestrat(Sequence seq, Protein* proteins, curandState* states) {        
    // curand_init(345678, threadIdx.x, 0, &states[threadIdx.x]);(
    // int cost = curanda(&states[threadIdx.x]);
    int cost = get_cost(*(proteins+threadIdx.x));
    for (int i = 0; i < 100; i++) {
		step(proteins[threadIdx.x], &states[threadIdx.x]);
		cost = get_cost(*(proteins+threadIdx.x));
    }
    printf("%d %d\n", threadIdx.x, cost);
    for (int offset = 32; offset > 0; offset /= 2) {
        auto other = __shfl_down_sync(0xFFFFFFFF, cost, offset);
        cost = cost < other ? cost : other;
    }
    printf("%d %d\n", threadIdx.x, cost);
}


void anneal_multistart_singlestrat(Sequence seq) {    
    Protein* proteins;	
    cudaMallocManaged(&proteins, sizeof(Protein) * 32);	
    curandState *dev_random;
    cudaMalloc((void**)&dev_random, 32*sizeof(curandState));
    for (int i = 0; i < 32; i++) proteins[i] = Protein::random(seq);
    __anneal_multistart_singlestrat<<<1,32>>>(seq, proteins, dev_random);
    cudaDeviceSynchronize();	
}

