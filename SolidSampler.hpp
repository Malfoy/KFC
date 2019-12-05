#ifndef __SOLID_SAMPLER__
#define __SOLID_SAMPLER__

#include "blocked.hpp"
#include <vector>

using kmer_t = size_t;

class SolidSampler {
    static constexpr size_t counter_bits = 2;
    static constexpr size_t counter_hashes = 3;

    static constexpr size_t saved_bloom_hashes = 10;
    static constexpr size_t saved_bloom_bits_per_kmer = saved_bloom_hashes;
    // Half of the memory is reserved for bits_per_kmer * m_kmer_max
    static constexpr size_t bits_per_kmer = saved_bloom_bits_per_kmer + CHAR_BIT * sizeof(kmer_t);

    using cbf_t = BlockedCMS<counter_bits, counter_hashes>;
    using saved_t = BlockedBloom<saved_bloom_hashes>;

  public:
    SolidSampler(size_t memory_size);

	void clean();
    bool insert(kmer_t kmer);
    std::vector<kmer_t>& get_kmers();

	friend std::ostream& operator<<(std::ostream& out, SolidSampler&);

    size_t get_nb_kmer_seen() const { return m_nb_inserted; }

  private:
    const size_t m_memory_cbf;
    const size_t m_kmer_max;
    const size_t m_memory_saved;
    size_t m_nb_inserted = 0;
    cbf_t m_cbf;
    saved_t m_saved;
    std::vector<kmer_t> kmers;
};

#endif
