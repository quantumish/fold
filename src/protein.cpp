#include <stdexcept>
#include <iostream>
#include <random>
#include "protein.hpp"

using namespace cpu;

Sequence::Sequence() {
    contents = nullptr;
	size = 0;
}

Sequence::Sequence(const std::string& str) {
    std::string raw_in = str;
    // std::transform(str.begin(), str.end(), raw_in.begin(), [](char c) {return std::toupper(c);});    
    size = str.length();
    for (size_t i = 0; i < str.length(); i++) {
        switch (str[i]) {
        case 'A': contents.push_back(Amino::Alanine); break;
        case 'R': contents.push_back(Amino::Arganine); break;
        case 'N': contents.push_back(Amino::Asparagine); break;
        case 'D': contents.push_back(Amino::AsparticAcid); break;
        case 'C': contents.push_back(Amino::Cysteine); break;
        case 'Q': contents.push_back(Amino::Glutamine); break;
        case 'E': contents.push_back(Amino::GlutamicAcid); break;
        case 'G': contents.push_back(Amino::Glycine); break;
        case 'H': contents.push_back(Amino::Histidine); break;
        case 'I': contents.push_back(Amino::Isoleucine); break;
        case 'L': contents.push_back(Amino::Leucine); break;
		case 'K': contents.push_back(Amino::Lysine); break;
        case 'M': contents.push_back(Amino::Methionine); break;
		case 'F': contents.push_back(Amino::Phenylalanine); break;
        case 'P': contents.push_back(Amino::Proline); break;
        case 'S': contents.push_back(Amino::Serine); break;
        case 'T': contents.push_back(Amino::Threonine); break;
        case 'W': contents.push_back(Amino::Tryptophan); break;
        case 'Y': contents.push_back(Amino::Tyrosine); break;
        case 'V': contents.push_back(Amino::Valine); break;
        default: throw std::domain_error("Invalid amino acid!");
        }
    }    
}

Protein::Protein(Sequence seq, Eigen::Vector3i* pos) :sequence{seq}, positions{pos} {}

// Initializes `Protein` in denatured state (a straight line across the Y axis)
Protein::Protein(Sequence seq) {
    sequence = seq;
	positions = new Eigen::Vector3i[seq.size];      
    for (int i = 0; i < seq.size; i++) {
        positions.push_back({0, i, 0};
    }
}

// Initializes 'Protein' in a random configuration. Useful for multistart methods.
Protein Protein::random(Sequence seq) {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> distr(0, 2);
    Eigen::Vector3i* pos = new Eigen::Vector3i[seq.size];
    pos[0] = {0,0,0};
    for (int i = 0; i < seq.size-1; i++) {
        bool clear = false;
        Eigen::Vector3i candidate;
        while (clear == false) {
            candidate = pos[i];			
            candidate[distr(gen)] += 1;
            // Ensure there are no overlapping candidates.
            clear = true;
            for (int j = 0; j < i; j++) {
                if (pos[j] == candidate) {
                    clear = false;
                }
            }
        }
        pos[i+1] = candidate;
    }
    return {seq, pos};    
}
