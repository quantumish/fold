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

#include "protein.hpp"

__device__ void step(Protein& protein, curandState* state) {    
    auto index = curand(state) % protein.sequence.size;
	auto pos = protein.positions[index];	
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
    for (int offset = 32; offset > 0; offset /= 2) {
        auto other = __shfl_down_sync(0xFFFFFFFF, cost, offset);
        cost = cost < other ? cost : other;
    }
}

void copy_protein(Protein& p, void* buffer) {	
	cudaMallocManaged(&p.positions, p.sequence.size*sizeof(Eigen::Vector3i));
	cudaMallocManaged(&p.sequence.contents, p.sequence.size*sizeof(Amino));
	cudaMemcpy(&p, buffer, sizeof(Protein));
}

void anneal_multistart_singlestrat(Sequence seq) {
	int rank, nprocs;
	MPI_Init(0, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	Protein* proteins = new Protein[32];	
    for (int i = 0; i < 32; i++) proteins[i] = Protein::random(seq);
    Protein* dev_proteins;	
    cudaMallocManaged(&dev_proteins, sizeof(Protein) * 32);	
	for (int i = 0; i < 32; i++) dev_proteins[i] = copy_protein(proteins[i]);
    curandState *dev_random;
    cudaMalloc((void**)&dev_random, 32*sizeof(curandState));
    __anneal_multistart_singlestrat<<<1,32>>>(seq, proteins, dev_random);
    while (true) { 
		if (MPI
	}
}

