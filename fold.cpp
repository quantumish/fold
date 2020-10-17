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
    void update();
};

Protein::Protein(std::string sequence)
{
    for (char c : sequence) {
        //if (polar.find(c) > 0) residues.push_back({c, true, {0,residues.size()}});
        //else if (nonpolar.find(c) > 0) residues.push_back({c, false, {0,residues.size()}});
        if (c = 'P') residues.push_back({c, true, {0,residues.size()}});
        else residues.push_back({c, false, {0,residues.size()}});
    }
    energy = 0;
};

Protein::update()
{
    
};

// Takes a index in a residue chain, determines what moves are valid for it
// Returns vector of bool where index 0 represents if end move is valid.
// index 1 represents corner flip, and index 2 represents crankshaft
std::vector<bool> check_moves(std::vector<Residue> residues, int index)
{
    std::vector<bool> valid;
    if (i == 0 || i == residues.size()) {
        if (check_residue({residue[i].coords[0]+1, residue[i].coords[1]}) &&
            check_residue({residue[i].coords[0]-1, residue[i].coords[1]}) &&
            check_residue({residue[i].coords[0], residue[i].coords[1]+1}) &&
            check_residue({residue[i].coords[0], residue[i].coords[1]-1})) {
            valid.push_back(false);
        }
        else valid.push_back(true)
    }
    
}

// Takes residue chain and a coordinate, checks if there is a residue there. O(n).
bool check_residue(std::vector<Residue> residues, Eigen::Vector2i target)
{
    for (int i = 0; i < residues.size(); i++) {
        if (residues[i].coords == target) return true;
    }
    return false;
}

// Takes residue chain and calculates energy score. O(n^2).
int energy(std::vector<Residue> residues)
{
    int energy = 0;
    for (int i = 0; i < residues.size(); i++) {
        if (residues[i].polar == true) continue;
        for (int j = 0; j < residues.size(); j++) {
            if (residues[j].polar == true || abs(i-j) == 1) continue;
            if (abs(residues[i].coords[0]-residues[j].coords[0]) == 1 ||
                abs(residues[i].coords[1]-residues[j].coords[1]) == 1) {
                energy-=1;
            }
        }
    }
    return energy/2;
}

int main()
{
    Protein protein ("HPPH");
    std::cout << energy(protein.residues) << "\n";
}
