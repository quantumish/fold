#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <string>
#include <cstdlib>

std::string polar = "QNHSTYC";
std::string nonpolar = "AILMFVPG";

struct Residue
{
    char id;
    bool polar;
    Eigen::Vector2i coords;
};

struct Protein
{
    std::vector<Residue> residues;
    int energy;
    Protein(std::string sequence);
    void fold();
};

Protein::Protein(std::string sequence)
{
    for (char c : sequence) {
        if (polar.find(c) > 0) residues.push_back({c, true, {0,residues.size()}});
        else if (nonpolar.find(c) > 0) residues.push_back({c, false, {0,residues.size()}});
    }
    energy = 0;
};

int energy(std::vector<Residue> protein)
{
    int energy = 0;
    for (int i = 0; i < protein.size(); i++) {
        if (protein[i].polar == true) continue;
        for (int j = 0; j < protein.size(); j++) {
            if (protein[j].polar == true || abs(i-j) == 1) continue;
            if (abs(protein[i].coords[0]-protein[j].coords[0]) == 1 ||
                abs(protein[i].coords[1]-protein[j].coords[1]) == 1) {
                energy-=1;
            }
        }
    }
    return energy/2;
}

int main()
{
    Protein protein ("ANNA");
    std::cout << energy(protein.residues) << "\n";
}
