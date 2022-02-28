#pragma once
#include <string>
#include <Eigen/Dense>

enum class Amino {
    Alanine,
    Arganine,
    Asparagine,
    AsparticAcid,
    Cysteine,
    Glutamine,
    GlutamicAcid,
    Glycine,
    Histidine,
    Isoleucine,
    Leucine,
    Lysine,
    Methionine,
    Phenylalanine,
    Proline,
    Serine,
    Threonine,
    Tryptophan,
    Tyrosine,
    Valine,
};

class Sequence {
	Sequence();
public:
	Amino* contents;
	size_t sz;
	explicit Sequence(const std::string& str);
};

class Protein {
	Protein(Sequence seq, Eigen::Vector3i* pos);
public:
    Amino* sequence;
	Eigen::Vector3i* positions;
	size_t sz;
	explicit Protein(Sequence seq);
	static Protein random(Sequence seq);
};
