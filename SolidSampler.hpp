#ifndef __SOLID_SAMPLER__
#define __SOLID_SAMPLER__

#include "blocked.hpp"
#include "CascadingBloomFilter.hpp"

class SolidSampler {
  public:
	SolidSampler(uint64_t memory_size);
	~SolidSampler();

	void clean();
	bool insert(uint8_t* kmer, std::size_t len);
	std::vector<uint64_t>& get_kmers();

	friend std::ostream& operator<<(std::ostream& out, SolidSampler&);

  private:
	uint64_t memory_cbf;
	BlockedCMS<> m_cbf;
	uint64_t m_nb_inserted;
	const uint64_t m_kmer_max;
	uint64_t m_nb_kmers_saved;
	std::vector<uint64_t> kmers;
	bool alive;
};

#endif
