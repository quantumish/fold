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
	auto index = curand(state) % protein.sz;
	auto pos = protein.positions[index];
}

__device__ bool check_residue(Protein protein, Eigen::Vector3i pos) {
	for (int i = 0; i < protein.sz; i++) {
		if (protein.positions[i] == pos) return false;
	}
	return true;
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
}

__global__ void __anneal_multistart_singlestrat(Sequence seq, Protein* proteins, int* min_cost, curandState* states) {
	int cost = get_cost(*(proteins+threadIdx.x));
	for (int i = 0; i < 100; i++) {
		step(proteins[threadIdx.x], &states[threadIdx.x]);
		cost = get_cost(*(proteins+threadIdx.x));
	}
	for (int offset = 32; offset > 0; offset /= 2) {
		auto other = __shfl_down_sync(0xFFFFFFFF, cost, offset);
		cost = cost > other ? cost : other;
	}	
	if (threadIdx.x == 0) *min_cost = cost;	
}

// void copy_protein(Protein& p, void* buffer) {
//	cudaMallocManaged(&p.positions, p.sequence.size*sizeof(Eigen::Vector3i));
//	cudaMallocManaged(&p.sequence, p.sequence.size*sizeof(Amino));
// }

void anneal_multistart_singlestrat(Sequence seq) {
	// int rank, nprocs;
	// MPI_Init(0, NULL);
	// MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	Protein* dev_proteins;
	cudaMallocManaged(&dev_proteins, sizeof(Protein) * 32);
	for (int i = 0; i < 32; i++) dev_proteins[i] = Protein::random(seq);
	curandState *dev_random;
	cudaMalloc(&dev_random, 32*sizeof(curandState));	
	int* cost;
	cudaMalloc(&cost, sizeof(int));
	__anneal_multistart_singlestrat<<<1,32>>>(seq, dev_proteins, cost, dev_random);
	cudaDeviceSynchronize();
	int* local_cost = new int;
	cudaMemcpy(cost, local_cost, 4, cudaMemcpyDeviceToHost);
	std::cout << *local_cost << "\n";
}
