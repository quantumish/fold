#include <Eigen/Dense>
#include <vector>
#include <array>
#include <iostream>
#include <optional>
#include <string>
#include <cstdlib>
#include <random>
#include <ctime>
#include <chrono>
#include <map>

#ifdef PYTHON
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
namespace py = pybind11;
#endif

constexpr double BOLTZ_CONST = 1.380649e-23;
std::map<char, int> amino_acid_map = {
    {'C', 0},  {'M', 1},  {'F', 2},  {'I', 3},  {'L', 4},  {'V', 5},  {'W', 6},
    {'Y', 7},  {'A', 8},  {'G', 9},  {'T', 10}, {'S', 11}, {'N', 12}, {'Q', 13},
    {'D', 14}, {'E', 15}, {'H', 16}, {'R', 17}, {'K', 18}, {'P', 19},
};

float interactions[20][20] = {
	{-3.477, -2.240, -2.424, -2.410, -2.343, -2.258, -2.080, -1.892, -1.700, -1.101,  -1.243, -1.306, -0.788, -0.835, -0.616, -0.179, -1.499, -0.771, -0.112, -1.196},
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
	{-1.196, -0.788, -1.076, -0.991, -0.771, -0.886, -1.278, -1.067, -0.374, -0.042,  -0.222, -0.199, -0.018, -0.035,  0.189,  0.257, -0.346, -0.023,  0.661,  0.129}
};

struct Residue
{
	int id;	
	Eigen::Vector3i backbone;
	Eigen::Vector3i side_chain;
};

/*
Takes residue chain `residues` and three-dimensional coordinates
and checks if there is a residue or sidechain at supplied coordinates.
If there is one there, returns the location of the residue in the protein sequence.
If there is no such residue, return empty std::nullopt. O(n) 
*/
static std::optional<size_t> check_for_entity(std::vector<Residue> residues, Eigen::Vector3i loc)
{	
	for (int i = 0; i < residues.size(); i++) {
		// FIXME Addition of sidechains makes logic very questionable.
		if (residues[i].backbone == loc || residues[i].backbone+residues[i].side_chain == loc) return i;
	}
    return std::nullopt;
}

// Overloaded `check_for_entity()` that adds a given 3d offset to
// given coordinates yet otherwise acts the same.
static std::optional<size_t> check_for_entity(std::vector<Residue> residues, Eigen::Vector3i loc, Eigen::Vector3i offset)
{
	loc += offset;
	for (int i = 0; i < residues.size(); i++) {
		if (residues[i].backbone == loc || residues[i].backbone+residues[i].side_chain == loc) return i;
	}
    return std::nullopt;
}

int exposure(std::vector<Residue> residues, size_t i)
{
	int n = 0;
	if (!check_for_entity(residues, residues[i].backbone, {-1, 0, 0}).has_value()) n++;
	if (!check_for_entity(residues, residues[i].backbone, {1, 0, 0}).has_value()) n++;
	if (!check_for_entity(residues, residues[i].backbone, {0, -1, 0}).has_value()) n++;
	if (!check_for_entity(residues, residues[i].backbone, {0, 1, 0}).has_value()) n++;
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
			if (((residues[i].side_chain+residues[i].backbone)-(residues[j].side_chain+residues[i].backbone)).norm() == 1 &&
				((residues[i].side_chain == residues[j].side_chain) || (residues[i].side_chain == -residues[j].side_chain))) {
				energy+=interactions[residues[i].id][residues[j].id];
			}
		}
	}
	return energy/2;
}

class Protein
{
	void find_end_moves(std::vector<std::function<void(void)>>& updates, size_t i);
	void find_corner_moves(std::vector<std::function<void(void)>>& updates, size_t i);
	void find_sidechain_moves(std::vector<std::function<void(void)>>& updates, size_t i);
public:
	std::vector<Residue> residues;
	float score;
	float temperature;
	Protein(std::string sequence, float temp, bool denatured);
	void update();
	int attempt_move(size_t i);
	int exposure(size_t i);
};

Protein::Protein(std::string sequence, float temp, bool denatured)
	:temperature(temp), score(0)
{
	for (int i = 0; i < sequence.size(); i++) {
		Eigen::Vector3i loc;
		// TODO Address the fact that this is a narrowing conversion from size_t -> int
		if (denatured) loc = {0,residues.size(),0};
		else if (i==0) loc = {0,0,0};
		else {			
			std::vector<Eigen::Vector3i> open_offsets;
			for (int j = 0; j < 3; j++) {
				Eigen::Vector3i offset = {0,0,0};
				for (int k : {1, -1}) {
					offset[j] = k;
					if (!check_for_entity(residues, residues[i-1].backbone, offset).has_value()) {
						open_offsets.push_back(residues[i-1].backbone+offset);
					}
				}
			}
			if (open_offsets.size() == 0) {
				throw std::runtime_error("Random protein formation failed: no possible locations for next residue.");
			}
			loc = open_offsets[rand() % open_offsets.size()];
		}
		residues.push_back({static_cast<int>(amino_acid_map[sequence[i]]), loc, loc+Eigen::Vector3i(1,0,0)});
	}
};

void Protein::update()
{
	std::vector<Residue> old = residues;
	auto now = std::chrono::high_resolution_clock::now();
	int code = attempt_move(rand() % residues.size());
	if (code == -1) return;
	float updated = energy(residues);
	float delta = updated - energy(old);
	float rndm = (float)rand()/RAND_MAX;
	if (rndm > sigmoidish(delta, temperature)) residues = old;
	else score = updated;
};

void Protein::find_end_moves(std::vector<std::function<void(void)>>& updates, size_t i)
{
	if (i == 0 || i == residues.size()-1) {
		int prev = 1;
		if (i == residues.size()-1) prev = -1;
		bool end = false;
		for (int j = 0; j < 3; j++) {
			for (int k : {1, -1}) {
				Eigen::Vector3i offset = {0,0,0};
				offset[j] = k;
				if (!check_for_entity(residues, residues[i+prev].backbone, offset).has_value()) {
					updates.emplace_back([this,i,prev,offset](){residues[i].backbone = residues[i+prev].backbone+offset;});
				}
			}
		}
	}
}


void Protein::find_corner_moves(std::vector<std::function<void(void)>>& updates, size_t i)
{
	// List of corner offsets for checking and generating fold offsets in a for loop (instead of 20 lines of if statements).
	// The amount of braces required for something like this is just absurd.
	// TODO: This is somewhat unacceptable.
	std::array<std::array<Eigen::Vector3i, 2>, 8> corners = {
		{{{{1,0,0},{0,-1,0}}}, {{{-1,0,0},{0,-1,0}}}, {{{-1,0,0},{0,1,0}}}, {{{1,0,0},{0,1,0}}},
		 {{{0,1,0},{0,0,-1}}}, {{{0,-1,0},{0,0,-1}}}, {{{0,-1,0},{0,0,1}}}, {{{0,1,0},{0,0,1}}}}
	};
	for (int j = 0; j < 8; j++) {
		// TODO/HACK: Review and simplify this if statement: it's likely buggy
		if (abs(static_cast<int>(i) - static_cast<int>(check_for_entity(residues, residues[i].backbone, corners[j][0]).value_or(NULL))) == 1 &&
			abs(static_cast<int>(i) - static_cast<int>(check_for_entity(residues, residues[i].backbone, corners[j][1]).value_or(NULL))) == 1) {			
			Eigen::Vector3i offset = corners[j][0]+corners[j][1];
			if (!check_for_entity(residues, residues[i].backbone, offset).has_value()) {
				updates.emplace_back([this,i,offset](){residues[i].backbone+=offset;});
				break;				
			}
		}
	}
}

void Protein::find_sidechain_moves(std::vector<std::function<void(void)>>& updates, size_t i)
{
	for (int j = 0; j < 3; j++) {
		for (int k : {1, -1}) {
			Eigen::Vector3i offset = {0,0,0};
			offset[j] = k;
			if (!check_for_entity(residues, residues[i].backbone, offset).has_value()) {
				updates.emplace_back([this,i,offset](){residues[i].side_chain = offset;});
			}
		}
	}
}


// Takes a index in the residue chain, tries a random move that is valid for it
int Protein::attempt_move(size_t i)
{
	std::vector<std::function<void(void)>> updates;
	find_end_moves(updates, i);
	find_corner_moves(updates, i);
    find_sidechain_moves(updates, i);
	if (updates.size() == 0) return -1;
	updates[rand() % updates.size()]();
	return 0;
}

int main()
{
	srand(time(0));
	Protein protein ("HPPHPH", 2, true);
	for (int i = 0; i < 30000; i++) {
		protein.update();
		// float cost = protein.energy;
	}
}

#ifdef PYTHON
PYBIND11_MODULE(fold, m) {
	m.doc() = "Protein folding.";
	py::class_<Residue>(m, "Residue")
		.def_readonly("id", &Residue::id)
		.def_readonly("sidechain", &Residue::side_chain)
		.def_readonly("backbone", &Residue::backbone);
	py::class_<Protein>(m, "Protein")
		.def(py::init<std::string, float, bool>())
		.def("update", &Protein::update)
		.def_readonly("residues", &Protein::residues)
		.def_readonly("score", &Protein::score);
	m.def("exposure", &exposure, py::arg("residues"), py::arg("i"));
}
#endif
