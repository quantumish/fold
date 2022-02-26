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

struct Sequence {
	Sequence();
public:
    Amino* contents;
    size_t size;	
    explicit Sequence(const std::string& str);
};


class Protein {
    Protein(Sequence seq, Eigen::Vector3i* pos);
public:
    Sequence sequence;
    Eigen::Vector3i* positions;
    explicit Protein(Sequence seq);
    static Protein random(Sequence seq);
};
