#ifndef __SOLID_SAMPLER__
#define __SOLID_SAMPLER__

#include "CascadingBloomFilter.hpp"

class SolidSampler {
  public:
	SolidSampler(uint64_t memory_size);

	void insert(uint8_t* kmer, std::size_t len);

	friend std::ostream& operator<<(std::ostream& out, SolidSampler&);

  private:
	CascadingBloomFilter* m_cbf;
	uint64_t m_nb_inserted;
	uint64_t m_kmer_max;
	std::vector<bool> saved;
	uint64_t m_nb_kmers_saved;
	std::vector<uint64_t> kmers;
};

#endif
