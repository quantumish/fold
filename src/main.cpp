#include <iostream>
#include "protein.hpp"
#include "fold.cuh"

int main() {
    Sequence seq ("HPHPHP");
	std::cout << seq.size << "\n";
	anneal_multistart_singlestrat(seq);
}
