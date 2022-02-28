#pragma once
#include <string>
#include <Eigen/Dense>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

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

namespace cpu {
	class Sequence {
		Sequence();
	public:
		thrust::host_vector<Amino> contents;    
		explicit Sequence(const std::string& str);
	};

	class Protein {
		Protein(Sequence seq, Eigen::Vector3i* pos);
	public:
		Sequence sequence;
		thrust::host_vector<Eigen::Vector3i> positions;
		explicit Protein(Sequence seq);
		static Protein random(Sequence seq);
	};
}

namespace gpu {
	struct Protein {
		thrust::device_vector<Amino> sequence;
		thrust::device_vector<Eigen::Vector3i> positions;
	};
}
