#include <stdexcept>
#include <iostream>
#include <random>
#include "protein.hpp"

Sequence::Sequence() {
    contents = nullptr;
	sz = 0;
}

Sequence::Sequence(const std::string& str) {
    std::string raw_in = str;
    // std::transform(str.begin(), str.end(), raw_in.begin(), [](char c) {return std::toupper(c);});    
    contents = new Amino[str.length()];
    sz = str.length();
    for (size_t i = 0; i < str.length(); i++) {
        switch (str[i]) {
        case 'A': contents[i] = Amino::Alanine; break;
        case 'R': contents[i] = Amino::Arganine; break;
        case 'N': contents[i] = Amino::Asparagine; break;
        case 'D': contents[i] = Amino::AsparticAcid; break;
        case 'C': contents[i] = Amino::Cysteine; break;
        case 'Q': contents[i] = Amino::Glutamine; break;
        case 'E': contents[i] = Amino::GlutamicAcid; break;
        case 'G': contents[i] = Amino::Glycine; break;
        case 'H': contents[i] = Amino::Histidine; break;
        case 'I': contents[i] = Amino::Isoleucine; break;
        case 'L': contents[i] = Amino::Leucine; break;
        case 'K': contents[i] = Amino::Lysine; break;
        case 'M': contents[i] = Amino::Methionine; break;
        case 'F': contents[i] = Amino::Phenylalanine; break;
        case 'P': contents[i] = Amino::Proline; break;
        case 'S': contents[i] = Amino::Serine; break;
        case 'T': contents[i] = Amino::Threonine; break;
        case 'W': contents[i] = Amino::Tryptophan; break;
        case 'Y': contents[i] = Amino::Tyrosine; break;
        case 'V': contents[i] = Amino::Valine; break;
        default: throw std::domain_error("Invalid amino acid!");
        }
    }    
}

Protein::Protein(Sequence seq, Eigen::Vector3i* pos) :sequence{seq.contents}, positions{pos}, sz{seq.sz} {}

// Initializes `Protein` in denatured state (a straight line across the Y axis)
Protein::Protein(Sequence seq)
	:sequence{seq.contents}, positions{new Eigen::Vector3i[seq.sz]}, sz{seq.sz}
{    
    for (size_t i = 0; i < seq.sz; i++) {
        positions[i] = {0, i, 0};
    }
}

// Initializes 'Protein' in a random configuration. Useful for multistart methods.
Protein Protein::random(Sequence seq) {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> distr(0, 2);
    Eigen::Vector3i* pos = new Eigen::Vector3i[seq.sz];
    pos[0] = {0,0,0};
    for (size_t i = 0; i < seq.sz-1; i++) {
        bool clear = false;
        Eigen::Vector3i candidate;
        while (clear == false) {
            candidate = pos[i];			
            candidate[distr(gen)] += 1;
            // Ensure there are no overlapping candidates.
            clear = true;
            for (size_t j = 0; j < i; j++) {
                if (pos[j] == candidate) {
                    clear = false;
                }
            }
        }
        pos[i+1] = candidate;
    }
    return {seq, pos};    
}
