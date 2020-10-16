#include <Eigen/Dense>
#include <cstdlib>

enum AminoAcid {A, R, N, D, C, E, Q, G, H, I, L, K, M, F, P, S, T, W, Y, V};

struct Residue
{
    AminoAcid id;
    bool polar;
    Eigen::Vector2i coords;
};

int energy(std::vector<Residue> protein)
{
    int energy = 0;
    for (int i = 0; i < protein.size(); i++) {
        if (protein[i].polar == true) continue;
        for (int j = 0; j < protein.size(); j++) {
            if (protein[i].polar == true) continue;
            if (abs(protein[i].coords[0]-protein[j].coords[0]) == 1 ||
                abs(protein[i].coords[1]-protein[j].coords[1]) == 1) {
                energy-=1;
            }
        }
    }
    return energy;
}

int main()
{
    std::vector<Residue> protein;
}
