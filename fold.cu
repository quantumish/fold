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

enum class Amino {
    Alanine,
    Arganine,
    Asparagine,
    AsparticAcid,
    Cysteine,
    Glutamine,
    GlutamicAcid,
    Glycine,
    Histidine,
    Isoleucine,
    Leucine,
    Lysine,
    Methionine,
    Phenylalanine,
    Proline,
    Serine,
    Threonine,
    Tryptophan,
    Tyrosine,
    Valine,
};

struct Sequence {
    Amino* contents;
    size_t size;
    Sequence(const std::string str);
};

Sequence::Sequence(const std::string str) {
    std::string raw_in = str;
    // std::transform(str.begin(), str.end(), raw_in.begin(), [](char c) {return std::toupper(c);});
    Sequence seq;
    cudaMallocManaged(&seq.contents, str.length());
    seq.size = str.length();
    for (int i = 0; i < str.length(); i++) {
        switch (str.get(i)) {
        case 'A': seq.contents[i] = Amino::Alanine; break;
        case 'R': seq.contents[i] = Amino::Arganine; break;
        case 'N': seq.contents[i] = Amino::Asparagine; break;
        case 'D': seq.contents[i] = Amino::AsparticAcid; break;
        case 'C': seq.contents[i] = Amino::Cysteine; break;
        case 'Q': seq.contents[i] = Amino::Glutamine; break;
        case 'E': seq.contents[i] = Amino::GlutamicAcid; break;
        case 'G': seq.contents[i] = Amino::Glycine; break;
        case 'H': seq.contents[i] = Amino::Histidine; break;
        case 'I': seq.contents[i] = Amino::Isoleucine; break;
        case 'L': seq.contents[i] = Amino::Leucine; break;
        case 'K': seq.contents[i] = Amino::Lysine; break;
        case 'M': seq.contents[i] = Amino::Methionine; break;
        case 'F': seq.contents[i] = Amino::Phenylalanine; break;
        case 'P': seq.contents[i] = Amino::Proline; break;
        case 'S': seq.contents[i] = Amino::Serine; break;
        case 'T': seq.contents[i] = Amino::Threonine; break;
        case 'W': seq.contents[i] = Amino::Tryptophan; break;
        case 'Y': seq.contents[i] = Amino::Tyrosine; break;
        case 'V': seq.contents[i] = Amino::Valine; break;
        default: throw std::domain_error("Invalid amino acid!");
        }
    }    
    return seq;
}

class Protein {
    Protein(Sequence seq, std::vector<Eigen::Vector3i> pos);
public:
    Sequence sequence;
    Eigen::Vector3i* positions;
    Protein(Sequence seq);
    static Protein random(Sequence seq);
};

Protein::Protein(Sequence seq, Eigen::Vector3i* pos) :sequence{seq}, positions{pos} {}

// Initializes `Protein` in denatured state (a straight line across the Y axis)
Protein::Protein(Sequence seq) {
    sequence = seq;
    cudaMallocManaged(&positions, seq.size);    
    for (int i = 0; i < seq.size; i++) {
        positions[i] = {0, i, 0};
    }
}

// Initializes 'Protein' in a random configuration. Useful for multistart methods.
Protein Protein::random(Sequence seq) {
    Eigen::Vector3i* pos;
    cudaMallocManaged(&pos, seq.size);
    pos[0] = {0,0,0};
    for (int i = 0; i < seq.size-1; i++) {
        bool clear = false;
        Eigen::Vector3i candidate;
        while (clear == false) {
            candidate = pos.back();
            candidate[rand() % 3] += 1; // cuRAND ify me
            // Ensure there are no overlapping candidates.
            clear = true;
            for (int j = 0; j < i; j++) {
                if (pos[j] == candidate) {
                    clear = false;
                }
            }
        }
        pos[i] = candidate;n
    }
    return {seq, pos};    
}

__device__ void step(Protein& protein, curandState* state) {    
    // rand() % 
}


__device__ bool check_residue(Protein protein, Eigen::Vector3i pos) {
    for (Eigen::Vector3i p : protein.positions) {
	if (p == pos) return false;
    }
    return true;
}

__device__ int get_cost(Protein protein) {
    int cost = 0;
    const Eigen::Vector3i unit[3] = {{0,0,1}, {0,1,0}, {1,0,0}};
    for (Eigen::Vector3i p : protein.positions) {
	for (int i = 0; i < 3; i++) {
	    if (check_residue(protein, p + unit[i])) cost++;
	    if (check_residue(protein, p - unit[i])) cost++;	    
	}
    }
}

__global__ void __anneal_multistart_singlestrat(Sequence seq, Protein* proteins, curandState* states) {        
    // curand_init(345678, threadIdx.x, 0, &states[threadIdx.x]);
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

int main() {
    auto seq = make_sequence("HPHPHP");
    anneal_multistart_singlestrat(seq);
}