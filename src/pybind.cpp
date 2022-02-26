#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
namespace py = pybind11;
#include "protein.hpp"

PYBIND11_MODULE(fold, m) {
  m.doc() = "Cheap low-resolution lattice protein folding simulation code.";
}
