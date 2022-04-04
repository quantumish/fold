#include <iostream>
#include "protein.hpp"
#include "fold.cuh"

int main() {
    Sequence seq ("HPHPHPHPHPHP");
    // while (true) {
    // 	Protein::random(seq);
    // }
    anneal_multistart_singlestrat(seq);
}

