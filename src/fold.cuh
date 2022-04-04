#pragma once
#include "protein.hpp"

std::tuple<std::vector<Protein>, int> anneal_multistart_singlestrat(Sequence seq);

