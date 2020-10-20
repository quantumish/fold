#include <Eigen/Dense>
#include <vector>
#include <array>
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

#define BOLTZ_CONST (1.380649 * pow(10, -23))

std::string polar = "QNHSTYC";
std::string nonpolar = "AILMFVPG";
enum Move {End, Corner, Crankshaft};

struct Residue
{
    char id;
    bool polar;
    Eigen::Vector2i coords;
};

// Takes residue chain and a coordinate, checks if there is a residue there. O(n).
int check_residue(std::vector<Residue> residues, Eigen::Vector2i target, Eigen::Vector2i offset)
{
    target += offset;
    //    std::cout << "Checking if there's a residue at " << target.transpose() << "...\n";
    for (int i = 0; i < residues.size(); i++) {
        if (residues[i].coords == target) return i;
    }
    //    std::cout << "Nope!" << "\n";
    return -5000;
}

int exposure(std::vector<Residue> residues, int i)
{
    int n = 0;
    if (check_residue(residues, residues[i].coords, {-1, 0}) < 0) n++;
    if (check_residue(residues, residues[i].coords, {1, 0}) < 0) n++;
    if (check_residue(residues, residues[i].coords, {0, -1}) < 0) n++;
    if (check_residue(residues, residues[i].coords, {0, 1}) < 0) n++;
    return n;
}

float sigmoidish(float x, float T)
{
    if (x>0) return 1;
    else return 1/(20+exp(x/T));
}

// Takes residue chain and calculates energy score. O(n^2).
float energy(std::vector<Residue> residues)
{
    float energy = 0;
    for (int i = 0; i < residues.size(); i++) {
        for ( int j = 0; j < residues.size(); j++) {
            if ((residues[j].polar && residues[i].polar) || abs(i-j) == 1 || i==j) continue;
            if ((abs(residues[i].coords[0]-residues[j].coords[0]) == 1 &&
                abs(residues[i].coords[1]-residues[j].coords[1]) == 0) ||
                (abs(residues[i].coords[0]-residues[j].coords[0]) == 0 &&
                abs(residues[i].coords[1]-residues[j].coords[1]) == 1)) {
                //if (residues[i].polar && !residues[j].polar) energy-=1.2;
                //else if (!residues[i].polar && residues[j].polar) energy-=1.2;
                if (!residues[i].polar && !residues[j].polar) energy-=1;
            }
        }
    }
    return energy/2;
}

struct Protein
{
    std::chrono::high_resolution_clock::time_point born;
    std::vector<Residue> residues;
    float score;
    float temperature;
    Protein(std::string sequence, float temp, bool denatured);
    void update();
    int attempt_move(int i);
    int exposure(int i);
};

Protein::Protein(std::string sequence, float temp, bool denatured)
    :temperature(temp)
{
    born = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < sequence.size(); i++) {
        Eigen::Vector2i loc;
        if (denatured) loc = {0,residues.size()};
        else if (i==0) loc = {0,0};
        else {
            std::vector<Eigen::Vector2i> open_offsets;
            std::array<Eigen::Vector2i, 4> offsets = {{{1,0},{-1,0},{0,1},{0,-1}}};
            for (Eigen::Vector2i offset : offsets) if (check_residue(residues, residues[i-1].coords, offset) < 0) open_offsets.push_back(residues[i-1].coords+offset);
            auto now = std::chrono::high_resolution_clock::now();
            srand(std::chrono::duration_cast<std::chrono::nanoseconds>(now - born).count());
            if (open_offsets.size() == 0) {
                std::cout << "Uh oh! We folded ourself into a corner." << "\n";
                exit(1);
            }
            loc = open_offsets[rand() % open_offsets.size()];
        }
        if (sequence[i] == 'P') residues.push_back({sequence[i], true, loc});
        else residues.push_back({sequence[i], false, loc});
    }
    score = 0;
};

void Protein::update()
{
    std::vector<Residue> old = residues;   
    auto now = std::chrono::high_resolution_clock::now();
    srand(std::chrono::duration_cast<std::chrono::nanoseconds>(now - born).count());
    int code = attempt_move(rand() % residues.size());
    if (code == -1) return;
    float updated = energy(residues);
    float delta = updated - energy(old);
    // now = std::chrono::high_resolution_clock::now();
    // srand(std::chrono::duration_cast<std::chrono::nanoseconds>(now - born).count());
    float rndm = (float)rand()/RAND_MAX;
    // std::cout << rndm << " vs exp((" << -delta << " + 0.1) /" << temperature << ") so " << exp((-delta + 0.1)/temperature) << "\n";
    if (rndm > sigmoidish(delta, temperature)) residues = old;
    else score = updated;
};

// Takes a index in the residue chain, tries a random move that is valid for it
int Protein::attempt_move(int i)
{
    std::vector<Move> valid;
    std::vector<Eigen::Vector2i> open_offsets;
    if (i == 0 || i == residues.size()-1) {
        std::array<Eigen::Vector2i, 4> offsets = {{{1,0},{-1,0},{0,1},{0,-1}}};
        int prev;
        if (i == 0) prev = 1;
        if (i == residues.size()-1) prev = -1;
        for (Eigen::Vector2i offset : offsets) {
            if (check_residue(residues, residues[i+prev].coords, offset) < 0) {
                valid.push_back(End);
                //std::cout << "There's nothing at " << (residues[i+prev].coords+offset).transpose() << "\n";
                open_offsets.push_back(residues[i+prev].coords - residues[i].coords + offset);
                //std::cout << "Move to " << (residues[i+prev].coords - residues[i].coords + offset).transpose() << "is valid!\n";
            }
        }
    }
    Eigen::Vector2i jump;
    // List of corner offsets for checking and generating fold offsets in a for loop (instead of 20 lines of if statements).
    // The amount of braces required for something like this is just absurd.
    std::array<std::array<Eigen::Vector2i, 2>, 4> corners = {{{{{1,0},{0,-1}}}, {{{-1,0},{0,-1}}}, {{{-1,0},{0,1}}}, {{{1,0},{0,1}}}}};
    // for (int j = 0; j < 4; j++) {
    //     std::cout << corners[j][0].transpose() << "and" << corners[j][1].transpose() << "\n";
    // }
    for (int j = 0; j < 4; j++) {
        if (abs(i - check_residue(residues, residues[i].coords, corners[j][0])) == 1 &&
            abs(i - check_residue(residues, residues[i].coords, corners[j][1])) == 1) {
            //std::cout << "As " << (residues[i].coords+corners[j][0]).transpose() << " (" << corners[j][0] << ") and " << (residues[i].coords+corners[j][1]).transpose() << " are next to us we're fine for corner move\n";
            jump = corners[j][0]+corners[j][1];
            //std::cout << "Our jump would be to " << (residues[i].coords+jump).transpose() << "\n";
            if (check_residue(residues, residues[i].coords, jump) < 0) {
                valid.push_back(Corner);
                //std::cout << "There's nothing there so we're good!" << "\n";
                break;
            }
        }
    }
    if (valid.size() == 0) return -1;
    auto now = std::chrono::high_resolution_clock::now();
    srand(std::chrono::duration_cast<std::chrono::nanoseconds>(now - born).count());
    Move action = valid[rand() % valid.size()];
    switch (action) {
    case Corner:
        residues[i].coords += jump;
        break;
    case End:
        now = std::chrono::high_resolution_clock::now();
        srand(std::chrono::duration_cast<std::chrono::nanoseconds>(now - born).count());
        residues[i].coords += open_offsets[rand() % open_offsets.size()];
        break;
    }
    return 0;
}

int main()
{
    Protein protein ("HPPHPH", 2, true);
    for (int i = 0; i < 300; i++) {
        protein.update();
        float cost = energy(protein.residues);
        std::cout << cost << "\n";
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
        .def(py::init<std::string, float, bool>())
        .def("update", &Protein::update)
        .def_readonly("residues", &Protein::residues)
        .def_readonly("score", &Protein::score);
    m.def("exposure", &exposure, py::arg("residues"), py::arg("i"));
}
#endif
