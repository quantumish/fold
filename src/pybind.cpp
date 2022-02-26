#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
namespace py = pybind11;
#include "protein.hpp"
#include "fold.cuh"

// class NiceSequence {
// };

// class NiceProtein {	
// };

PYBIND11_MODULE(fold, m) {
	m.doc() = "Cheap low-resolution lattice protein folding simulation code.";
	py::class_<Sequence>(m, "Sequence")
		.def(py::init<const std::string&>());
	py::class_<Protein>(m, "RawProtein")
		.def(py::init<Sequence>())
		.def_static("random", &Protein::random);	
	m.def("anneal_multistart_singlestrat", &anneal_multistart_singlestrat);
}
