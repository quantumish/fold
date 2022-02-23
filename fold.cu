#include <exception>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <random>
#include <Eigen/Dense>

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

using Sequence = std::vector<Amino>;
Sequence make_sequence(const std::string str) {
    std::string raw_in;
    std::transform(str.begin(), str.end(), raw_in.begin(), [](char c) {return std::toupper(c);});
    Sequence seq;
    for (char c : raw_in) {
        switch (c) {
        case 'A': seq.push_back(Amino::Alanine); break;
        case 'R': seq.push_back(Amino::Arganine); break;
        case 'N': seq.push_back(Amino::Asparagine); break;
        case 'D': seq.push_back(Amino::AsparticAcid); break;
        case 'C': seq.push_back(Amino::Cysteine); break;
        case 'Q': seq.push_back(Amino::Glutamine); break;
        case 'E': seq.push_back(Amino::GlutamicAcid); break;
        case 'G': seq.push_back(Amino::Glycine); break;
        case 'H': seq.push_back(Amino::Histidine); break;
        case 'I': seq.push_back(Amino::Isoleucine); break;
        case 'L': seq.push_back(Amino::Leucine); break;
        case 'K': seq.push_back(Amino::Lysine); break;
        case 'M': seq.push_back(Amino::Methionine); break;
        case 'F': seq.push_back(Amino::Phenylalanine); break;
        case 'P': seq.push_back(Amino::Proline); break;
        case 'S': seq.push_back(Amino::Serine); break;
        case 'T': seq.push_back(Amino::Threonine); break;
        case 'W': seq.push_back(Amino::Tryptophan); break;
        case 'Y': seq.push_back(Amino::Tyrosine); break;
        case 'V': seq.push_back(Amino::Valine); break;
        default: throw std::domain_error("Invalid amino acid!");
        }
    }   
}

class Protein {
    Protein(Sequence seq, std::vector<Eigen::Vector3i> pos);
    Sequence sequence;
    std::vector<Eigen::Vector3i> positions;
public:
    Protein(Sequence seq);
    static Protein random(Sequence seq);
};

Protein::Protein(Sequence seq, std::vector<Eigen::Vector3i> pos) :sequence{seq}, positions{pos} {}

// Initializes `Protein` in denatured state (a straight line across the Y axis)
Protein::Protein(Sequence seq) {
    sequence = seq;
    for (int i = 0; i < seq.size(); i++) {
        positions.push_back({0, i, 0});
    }
}

// Initializes 'Protein' in a random configuration. Useful for multistart methods.
Protein Protein::random(Sequence seq) {
    std::vector<Eigen::Vector3i> pos;
    pos.push_back({0,0,0});
    for (int i = 0; i < seq.size()-1; i++) {
        bool clear = false;
        Eigen::Vector3i candidate;
        while (clear == false) {
            candidate = pos.back();
            candidate[rand() % 3] += 1;
            // Ensure there are no overlapping candidates.
            clear = true;
            for (Eigen::Vector3i j : pos) {
                if (j == candidate) {
                    clear = false;
                }
            }
        }
        pos.push_back(candidate);
    }
    return {seq, pos};    
}

__global__ void __anneal_multistart_singlestrat(Sequence seq) {
   
}


void anneal_multistart_singlestrat(Sequence seq) {
    
}