#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <string>
#include <cstdlib>
#include <random>
#include <ctime>
#include <chrono>

#ifdef PYTHON
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
namespace py = pybind11;
#endif

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
float energy(std::vector<Residue> residues)
{
    float energy = 0;
    for (int i = 0; i < residues.size(); i++) {
        if (residues[i].polar == true) continue;
        for (int j = 0; j < residues.size(); j++) {
            if (residues[j].polar == true || abs(i-j) == 1 || i==j) continue;
            //std::cout << i << "(" << residues[i].coords.transpose() << ") vs " << j << "(" << residues[j].coords.transpose() << ") is " << abs((residues[i].coords - residues[j].coords).norm()) << "\n";
            energy -= 1.0/abs((residues[i].coords - residues[j].coords).norm());
        }
    }
    return energy/2;
}

struct Protein
{
    std::chrono::high_resolution_clock::time_point born = std::chrono::high_resolution_clock::now();
    std::vector<Residue> residues;
    float score;
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
    float updated = energy(residues);
    if (updated > energy(old)) residues = old;
    else score = updated;
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
    if (i == 0 || i == residues.size()-1) {
        int prev;
        if (i == 0) prev = 1;
        if (i == residues.size()-1) prev = -1;
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
    bool top_left = (abs(i - check_residue(residues, residues[i].coords, 1, 0)) == 1 &&
                       abs(i - check_residue(residues, residues[i].coords, 0, -1)) == 1);
    bool top_right = (abs(i - check_residue(residues, residues[i].coords, -1, 0)) == 1 &&
                       abs(i - check_residue(residues, residues[i].coords, 0, -1)) == 1);
    bool bottom_right = (abs(i - check_residue(residues, residues[i].coords, -1, 0)) == 1 &&
                          abs(i - check_residue(residues, residues[i].coords, 0, 1)) == 1);
    bool bottom_left = (abs(i - check_residue(residues, residues[i].coords, 1, 0)) == 1 &&
                          abs(i - check_residue(residues, residues[i].coords, 0, 1)) == 1);
    //if (i == 1) std::cout << bottom_left << "\n";
    Eigen::Vector2i jump;
    if (bottom_right) {
        jump = {-1, 1};
        if (check_residue(residues, residues[i].coords, -1, 1) < 0) valid.push_back(Corner);
    } else if (top_left) {
        jump = {1, -1};
        if (check_residue(residues, residues[i].coords, 1, -1) < 0) valid.push_back(Corner);
    } else if (top_right) {
        jump = {-1, -1};
        if (check_residue(residues, residues[i].coords, -1, -1) < 0) valid.push_back(Corner);
    } else if (bottom_left) {
        jump = {1, 1};
        //        std::cout << "hi" << "\n";
        //        std::cout << check_residue(residues, residues[i].coords, 1, 1) << "\n";
        if (check_residue(residues, residues[i].coords, 1, 1) < 0) {
            //  std::cout << "yay" << "\n";
            valid.push_back(Corner);
        }
    }
    if (valid.size() == 0) return;
    auto now = std::chrono::high_resolution_clock::now();
    srand(std::chrono::duration_cast<std::chrono::nanoseconds>(now - born).count());
    Move action = valid[rand() % valid.size()];
    switch (action) {
    case Corner:
        //        std::cout << residues[i].coords.transpose() << " " << jump.transpose() << "\n";
        residues[i].coords += jump;
        //        std::cout << residues[i].coords.transpose() << "\n";
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
    Protein protein ("HPPHPH");
    for (int i = 0; i < 300; i++) {
        protein.update();
        float cost = energy(protein.residues);
        std::cout << cost << "\n";
        // if (cost <= -2) {
        //     for (int j = 0; j < protein.residues.size(); j++) {
        //         if (protein.residues[j].polar == false) std::cout << "(H)" << "\n";
        //         std::cout << protein.residues[j].coords << "\n\n";
        //     }
        //     return 0;
        // }
    }
    for (int j = 0; j < protein.residues.size(); j++) {
        if (protein.residues[j].polar == false) std::cout << "(H)" << "\n";
        std::cout << protein.residues[j].coords << "\n\n";
    }
}

#ifdef PYTHON
PYBIND11_MODULE(fold, m) {
    m.doc() = "Protein folding.";
    py::class_<Residue>(m, "Residue")
        .def_readonly("id", &Residue::id)
        .def_readonly("coords", &Residue::coords)
        .def_readonly("polar", &Residue::polar);
    py::class_<Protein>(m, "Protein")
        .def(py::init<std::string>())
        .def("update", &Protein::update)
        .def_readonly("residues", &Protein::residues)
        .def_readonly("score", &Protein::score);
}
#endif
