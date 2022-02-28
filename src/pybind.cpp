#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
namespace py = pybind11;
#include "protein.hpp"
#include "fold.cuh"

auto get_residues(Protein p) {
	std::vector<std::tuple<Amino, Eigen::Vector3i>> out;
	for (int i = 0; i < p.sz; i++) {
		out.push_back({p.sequence[i], p.positions[i]});
	}
	return out;
}

PYBIND11_MODULE(fold, m) {
	m.doc() = "Cheap low-resolution lattice protein folding simulation code.";
	py::enum_<Amino>(m, "Amino")
		.value("Alanine", Amino::Alanine)
		.value("Arganine", Amino::Arganine)
		.value("Asparagine", Amino::Asparagine)
		.value("AsparticAcid", Amino::AsparticAcid)
		.value("Cysteine", Amino::Cysteine)
		.value("Glutamine", Amino::Glutamine)
		.value("GlutamicAcid", Amino::GlutamicAcid)
		.value("Glycine", Amino::Glycine)
		.value("Histidine", Amino::Histidine)
		.value("Isoleucine", Amino::Isoleucine)
		.value("Leucine", Amino::Leucine)
		.value("Lysine", Amino::Lysine)
		.value("Methionine", Amino::Methionine)
		.value("Phenylalanine", Amino::Phenylalanine)
		.value("Proline", Amino::Proline)
		.value("Serine", Amino::Serine)
		.value("Threonine", Amino::Threonine)
		.value("Tryptophan", Amino::Tryptophan)
		.value("Tyrosine", Amino::Tyrosine)
		.value("Valine", Amino::Valine)
		.export_values();
	py::class_<Sequence>(m, "Sequence")
		.def(py::init<const std::string&>());
	py::class_<Protein>(m, "RawProtein")
		.def(py::init<Sequence>())
		.def_static("random", &Protein::random);	
	m.def("anneal_multistart_singlestrat", &anneal_multistart_singlestrat);
	m.def("get_residues", &get_residues);
}
