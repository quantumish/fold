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
std::string amino = "CMFILVWYAGTSNQDEHRKP";
float interactions[20][20] =
    {{-3.477, -2.240, -2.424, -2.410, -2.343, -2.258, -2.080, -1.892, -1.700, -1.101,  -1.243, -1.306, -0.788, -0.835, -0.616, -0.179, -1.499, -0.771, -0.112, -1.196},
     {-2.240, -1.901, -2.304, -2.286, -2.208, -2.079, -2.090, -1.834, -1.517, -0.897,  -0.999, -0.893, -0.658, -0.720, -0.409, -0.209, -1.252, -0.611, -0.146, -0.788}, 
     {-2.424, -2.304, -2.467, -2.530, -2.491, -2.391, -2.286, -1.963, -1.750, -1.034,  -1.237, -1.178, -0.790, -0.807, -0.482, -0.419, -1.330, -0.805, -0.270, -1.076}, 
     {-2.410, -2.286, -2.530, -2.691, -2.647, -2.568, -2.303, -1.998, -1.872, -0.885,  -1.360, -1.037, -0.669, -0.778, -0.402, -0.439, -1.234, -0.854, -0.253, -0.991}, 
     {-2.343, -2.208, -2.491, -2.647, -2.501, -2.447, -2.222, -1.919, -1.728, -0.767,  -1.202, -0.959, -0.524, -0.729, -0.291, -0.366, -1.176, -0.758, -0.222, -0.771}, 
     {-2.258, -2.079, -2.391, -2.568, -2.447, -2.385, -2.097, -1.790, -1.731, -0.756,  -1.240, -0.933, -0.673, -0.642, -0.298, -0.335, -1.118, -0.664, -0.200, -0.886}, 
     {-2.080, -2.090, -2.286, -2.303, -2.222, -2.097, -1.867, -1.834, -1.565, -1.142,  -1.077, -1.145, -0.884, -0.997, -0.613, -0.624, -1.383, -0.912, -0.391, -1.278}, 
     {-1.892, -1.834, -1.963, -1.998, -1.919, -1.790, -1.834, -1.335, -1.318, -0.818,  -0.892, -0.859, -0.670, -0.687, -0.631, -0.453, -1.222, -0.745, -0.349, -1.067}, 
     {-1.700, -1.517, -1.750, -1.872, -1.728, -1.731, -1.565, -1.318, -1.119, -0.290,  -0.717, -0.607, -0.371, -0.323, -0.235, -0.039, -0.646, -0.327,  0.196, -0.374}, 
     {-1.101, -0.897, -1.034, -0.885, -0.767, -0.756, -1.142, -0.818, -0.290,  0.219,  -0.311, -0.261, -0.230,  0.033, -0.097,  0.443, -0.325, -0.050,  0.589, -0.042}, 
     {-1.243, -0.999, -1.237, -1.360, -1.202, -1.240, -1.077, -0.892, -0.717, -0.311,  -0.617, -0.548, -0.463, -0.342, -0.382, -0.192, -0.720, -0.247,  0.155, -0.222}, 
     {-1.306, -0.893, -1.178, -1.037, -0.959, -0.933, -1.145, -0.859, -0.607, -0.261,  -0.548, -0.519, -0.423, -0.260, -0.521, -0.161, -0.639, -0.264,  0.223, -0.199}, 
     {-0.788, -0.658, -0.790, -0.669, -0.524, -0.673, -0.884, -0.670, -0.371, -0.230,  -0.463, -0.423, -0.367, -0.253, -0.344,  0.160, -0.455, -0.114,  0.271, -0.018}, 
     {-0.835, -0.720, -0.807, -0.778, -0.729, -0.642, -0.997, -0.687, -0.323,  0.033,  -0.342, -0.260, -0.253,  0.054,  0.022,  0.179, -0.290, -0.042,  0.334, -0.035}, 
     {-0.616, -0.409, -0.482, -0.402, -0.291, -0.298, -0.613, -0.631, -0.235, -0.097,  -0.382, -0.521, -0.344,  0.022,  0.179,  0.634, -0.664, -0.584, -0.176,  0.189}, 
     {-0.179, -0.209, -0.419, -0.439, -0.366, -0.335, -0.624, -0.453, -0.039,  0.443,  -0.192, -0.161,  0.160,  0.179,  0.634,  0.933, -0.324, -0.374, -0.057,  0.257}, 
     {-1.499, -1.252, -1.330, -1.234, -1.176, -1.118, -1.383, -1.222, -0.646, -0.325,  -0.720, -0.639, -0.455, -0.290, -0.664, -0.324, -1.078, -0.307,  0.388, -0.346}, 
     {-0.771, -0.611, -0.805, -0.854, -0.758, -0.664, -0.912, -0.745, -0.327, -0.050,  -0.247, -0.264, -0.114, -0.042, -0.584, -0.374, -0.307,  0.200,  0.815, -0.023}, 
     {-0.112, -0.146, -0.270, -0.253, -0.222, -0.200, -0.391, -0.349,  0.196,  0.589,   0.155,  0.223,  0.271,  0.334, -0.176, -0.057,  0.388,  0.815,  1.339,  0.661}, 
     {-1.196, -0.788, -1.076, -0.991, -0.771, -0.886, -1.278, -1.067, -0.374, -0.042,  -0.222, -0.199, -0.018, -0.035,  0.189,  0.257, -0.346, -0.023,  0.661,  0.129}};

struct Residue
{
    int id;
    bool polar;
    Eigen::Vector3i coords;
};

// Takes residue chain and a coordinate, checks if there is a residue there. O(n).
int check_residue(std::vector<Residue> residues, Eigen::Vector3i target, Eigen::Vector3i offset)
{
    target += offset;
    for (int i = 0; i < residues.size(); i++) {
        if (residues[i].coords == target) return i;
    }
    return -2;
}

int exposure(std::vector<Residue> residues, int i)
{
    int n = 0;
    if (check_residue(residues, residues[i].coords, {-1, 0, 0}) < 0) n++;
    if (check_residue(residues, residues[i].coords, {1, 0, 0}) < 0) n++;
    if (check_residue(residues, residues[i].coords, {0, -1, 0}) < 0) n++;
    if (check_residue(residues, residues[i].coords, {0, 1, 0}) < 0) n++;
    return n;
}

float sigmoidish(float x, float T)
{
    if (x<0) return 1;
    else return 1/(100+exp(x/T));
}

// Takes residue chain and calculates energy score. O(n^2).
float energy(std::vector<Residue> residues)
{
    float energy = 0;
    for (int i = 0; i < residues.size(); i++) {
        for ( int j = 0; j < residues.size(); j++) {
            if (abs(i-j) == 1 || i==j) continue;
            if ((residues[i].coords-residues[j].coords).norm() == 1) {
                energy+=interactions[residues[i].id][residues[j].id];
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
        Eigen::Vector3i loc;
        bool polarity = false;
        if (polar.find(sequence[i]) != std::string::npos) polarity = true;
        if (denatured) loc = {0,residues.size(),0};
        else if (i==0) loc = {0,0,0};
        else {
            std::vector<Eigen::Vector3i> open_offsets;
            for (int j = 0; j < 3; j++) {
                Eigen::Vector3i offset = {0,0,0};
                offset[j] = 1;
                if (check_residue(residues, residues[i-1].coords, offset) < 0) open_offsets.push_back(residues[i-1].coords+offset);
                offset[j] = -1;
                if (check_residue(residues, residues[i-1].coords, offset) < 0) open_offsets.push_back(residues[i-1].coords+offset);                
            }
            auto now = std::chrono::high_resolution_clock::now();
            srand(std::chrono::duration_cast<std::chrono::nanoseconds>(now - born).count());
            if (open_offsets.size() == 0) {
                std::cout << "Uh oh! We folded ourself into a corner." << "\n";
                exit(1);
            }
            loc = open_offsets[rand() % open_offsets.size()];
        }
        residues.push_back({static_cast<int>(amino.find(sequence[i])), polarity, loc});
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
    float rndm = (float)rand()/RAND_MAX;
    if (rndm > sigmoidish(delta, temperature)) residues = old;
    else score = updated;
};

// Takes a index in the residue chain, tries a random move that is valid for it
int Protein::attempt_move(int i)
{
    std::vector<Move> valid;
    std::vector<Eigen::Vector3i> open_offsets;
    if (i == 0 || i == residues.size()-1) {
        int prev = 1;
        if (i == residues.size()-1) prev = -1;
        bool end = false;
        for (int j = 0; j < 3; j++) {
            for (int k : {1, -1}) {
                Eigen::Vector3i offset = {0,0,0};
                offset[j] = k;
                if (check_residue(residues, residues[i+prev].coords, offset) < 0) {
                    end = true;
                    open_offsets.push_back(residues[i+prev].coords+offset);
                }
            }
        }
        if (end == true) valid.push_back(End);
    }
    Eigen::Vector3i jump;
    // List of corner offsets for checking and generating fold offsets in a for loop (instead of 20 lines of if statements).
    // The amount of braces required for something like this is just absurd.
    std::array<std::array<Eigen::Vector3i, 2>, 8> corners = {{{{{1,0,0},{0,-1,0}}}, {{{-1,0,0},{0,-1,0}}}, {{{-1,0,0},{0,1,0}}}, {{{1,0,0},{0,1,0}}},
                                                            {{{0,1,0},{0,0,-1}}}, {{{0,-1,0},{0,0,-1}}}, {{{0,-1,0},{0,0,1}}}, {{{0,1,0},{0,0,1}}}}};
    for (int j = 0; j < 8; j++) {
        if (abs(i - check_residue(residues, residues[i].coords, corners[j][0])) == 1 &&
            abs(i - check_residue(residues, residues[i].coords, corners[j][1])) == 1) {
            jump = corners[j][0]+corners[j][1];
            if (check_residue(residues, residues[i].coords, jump) < 0) {
                valid.push_back(Corner);
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
        //std::cout << "YAY CORNER [RES " << i << " MOVE BY " << jump.transpose() << "]\n";
        residues[i].coords += jump;
        break;
    case End:
        now = std::chrono::high_resolution_clock::now();
        srand(std::chrono::duration_cast<std::chrono::nanoseconds>(now - born).count());
        residues[i].coords = open_offsets[rand() % open_offsets.size()];
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
