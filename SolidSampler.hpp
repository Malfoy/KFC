#ifndef __SOLID_SAMPLER__
#define __SOLID_SAMPLER__

#include "CascadingBloomFilter.hpp"

class SolidSampler {
  public:
	SolidSampler(uint64_t memory_size);
	~SolidSampler();

	void clean();
	void insert(uint8_t* kmer, std::size_t len);
	std::vector<uint64_t>& get_kmers();

	friend std::ostream& operator<<(std::ostream& out, SolidSampler&);

  private:
	CascadingBloomFilter m_cbf;
	uint64_t m_nb_inserted;
	const uint64_t m_kmer_max;
	std::vector<bool> saved;
	uint64_t m_nb_kmers_saved;
	std::vector<uint64_t> kmers;
	bool alive;
	// std::vector<uint64_t> kmers;
};

#endif
