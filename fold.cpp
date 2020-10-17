#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <string>
#include <cstdlib>
#include <random>
#include <ctime>
#include <chrono>

std::string polar = "QNHSTYC";
std::string nonpolar = "AILMFVPG";
enum Move {End, Corner, Crankshaft};

struct Residue
{
    char id;
    bool polar;
    Eigen::Vector2i coords;
};

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

struct Protein
{
    std::chrono::high_resolution_clock::time_point born = std::chrono::high_resolution_clock::now();
    std::vector<Residue> residues;
    Protein(std::string sequence);
    void update();
    void attempt_move(int i);
};

Protein::Protein(std::string sequence)
{
    for (char c : sequence) {
        //if (polar.find(c) > 0) residues.push_back({c, true, {0,residues.size()}});
        //else if (nonpolar.find(c) > 0) residues.push_back({c, false, {0,residues.size()}});
        if (c == 'P') residues.push_back({c, true, {0,residues.size()}});
        else residues.push_back({c, false, {0,residues.size()}});
    }
};

void Protein::update()
{
    std::vector<Residue> old = residues;   
    auto now = std::chrono::high_resolution_clock::now();
    srand(std::chrono::duration_cast<std::chrono::nanoseconds>(now - born).count());
    attempt_move(rand() % residues.size());
    if (energy(residues) >= energy(old)) residues = old;
};

// Takes residue chain and a coordinate, checks if there is a residue there. O(n).
int check_residue(std::vector<Residue> residues, Eigen::Vector2i target, int x_offset = 0, int y_offset = 0)
{
    target[0] += x_offset;
    target[1] += y_offset;
    for (int i = 0; i < residues.size(); i++) {
        if (residues[i].coords == target) return i;
    }
    return -1;
}

// Takes a index in a residue chain, tries a random move that is valid for it
void Protein::attempt_move(int i)
{
    std::vector<Move> valid;
    std::vector<Eigen::Vector2i> open_offsets;
    std::cout << "i = " << i << "\n";
    if (i == 0 || i == residues.size()) {
        int prev;
        if (i == 0) prev = 1;
        if (i == residues.size()) prev = -1;
        if (check_residue(residues, residues[i+prev].coords, 1, 0) < 0) {
            valid.push_back(End);
            open_offsets.push_back({residues[i+prev].coords[0] - residues[i].coords[0] + 1, residues[i+prev].coords[1] - residues[i].coords[1] - 0});
        } else if (check_residue(residues, residues[i+prev].coords, -1, 0) < 0) {
            valid.push_back(End);
            open_offsets.push_back({residues[i+prev].coords[0] - residues[i].coords[0] - 1, residues[i+prev].coords[1] - residues[i].coords[1] - 0});
        } else if (check_residue(residues, residues[i+prev].coords, 0, 1) < 0) {
            valid.push_back(End);
            open_offsets.push_back({residues[i+prev].coords[0] - residues[i].coords[0] + 0, residues[i+prev].coords[1] - residues[i].coords[1] + 1});
        } else if (check_residue(residues, residues[i+prev].coords, 0, -1) < 0) {
            valid.push_back(End);
            open_offsets.push_back({residues[i+prev].coords[0] - residues[i].coords[0] + 0, residues[i+prev].coords[1] - residues[i].coords[1] - 1});
        }
    }
    bool top_corner = (abs(i - check_residue(residues, residues[i].coords, 1, 0)) == 1 &&
                       abs(i - check_residue(residues, residues[i].coords, 0, -1)) == 1);
    bool bottom_corner = (abs(i - check_residue(residues, residues[i].coords, -1, 0)) == 1 &&
                          abs(i - check_residue(residues, residues[i].coords, 0, 1)) == 1);
    Eigen::Vector2i jump;
    if (bottom_corner) {
        jump = {-1, 1};
        if (!check_residue(residues, residues[i].coords, -1, 1)) valid.push_back(Corner);
    } else if (top_corner) {
        jump = {1, -1};
        if (!check_residue(residues, residues[i].coords, 1, -1)) valid.push_back(Corner);
    }
    if (valid.size() == 0) return;
    auto now = std::chrono::high_resolution_clock::now();
    srand(std::chrono::duration_cast<std::chrono::nanoseconds>(now - born).count());
    Move action = valid[rand() % valid.size()];
    switch (action) {
    case Corner:
        residues[i].coords += jump;
        break;
    case End:
        auto now = std::chrono::high_resolution_clock::now();
        srand(std::chrono::duration_cast<std::chrono::nanoseconds>(now - born).count());
        residues[i].coords += open_offsets[rand() % open_offsets.size()];
        break;
    }
}

int main()
{
    Protein protein ("HPPH");
    for (int i = 0; i < 5; i++) {
        protein.update();
        int cost = energy(protein.residues);
        std::cout << "Energy: " << cost << "\n";
        for (int j = 0; j < protein.residues.size(); j++) {
            std::cout << protein.residues[j].coords << "\n\n";
        }
    }
}
