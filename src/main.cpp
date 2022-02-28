#include <iostream>
#include "protein.hpp"
#include "fold.cuh"

int main() {
	cpu::Sequence seq ("HPHPHP");
	anneal_multistart_singlestrat(seq);
}
