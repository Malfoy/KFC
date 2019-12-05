
#include <cmath>

#include <iostream>
#include <vector>
#include <algorithm>

#include "SolidSampler.hpp"
#include "Hash.hpp"
#include "blocked.hpp"
#include <gatbl/kmer.hpp>

uint8_t ilog2(size_t x) {
    return static_cast<uint8_t>(sizeof(size_t)) * CHAR_BIT - __builtin_clzll(x) - 1;
}

using namespace std;

/**
 * Construct the SolidSampler with maximum `memory_size` bytes for the cascading bloom filter.
 * For now, the solid kmers are not limited. So the SolidSampler can become huge.
 * @param cbf_memory_size Number of bytes to use in the CBF.
 * TODO : Allow a collision rate in the list => Compute the number of bits per kmer on the fly
 * and compute also the optimal number of hash functions.
 */
SolidSampler::SolidSampler(uint64_t memory_size)
  : m_memory_cbf(memory_size / 2)
  , m_kmer_max((CHAR_BIT * (memory_size - m_memory_cbf)) / bits_per_kmer)
  , m_memory_saved((this->m_kmer_max * saved_bloom_bits_per_kmer) / CHAR_BIT)
  , m_cbf(ilog2(m_memory_cbf / cbf_t::block_bytes))
  , m_saved(ilog2(m_memory_saved / saved_t::block_bytes))
  , kmers() {
    this->kmers.reserve(m_kmer_max);
}

/**
 * Insert a kmer into the Cascading Bloom Filter data structure.
 * If the kmer complete its path through the CBF, it's added into the frequent kmer vector.
 * @param kmer: The kmer wrapped into a byte array.
 * @param len: The byte array size.
 * @returns true when the number of saved kmer exceed the limit
 */
bool SolidSampler::insert(kmer_t kmer) {
    __uint128_t hash;
    MurmurHash3_x64_128(&kmer, sizeof(kmer_t), 0, &hash);

	this->m_nb_inserted++;

    if (!m_saved.possiblyContains(hash)) {
        auto old_count = m_cbf.add(hash);
        if (old_count == BlockedCMS<>::max_count) {
            // assert(std::find(kmers.begin(), kmers.end(), kmer) == kmers.end(), "Ouate le phoque");
            this->kmers.push_back(kmer);
            assume(m_saved.add(hash) == false, "kmer should not be present in the deduplicating bloom");
            return this->kmers.size() >= m_kmer_max;
        }
    }
	return false;
}

/**
 * Return the frequent kmers into a vector.
 */
vector<uint64_t>& SolidSampler::get_kmers() {
	return this->kmers;
}

ostream& operator<<(ostream& out, SolidSampler& sampler) {
    auto ncounters = sampler.m_cbf.size();
    auto saved_bits = sampler.m_saved.size();
    auto psolid = sampler.kmers.size();

    std::cout << "Inserted: " << sampler.m_nb_inserted << std::endl;

    out << "Counters:\n" << sampler.m_cbf;
    out << "\tmemory: " << double(ncounters * SolidSampler::counter_bits) / (size_t(CHAR_BIT) << 20) << "M\n";

    out << "Deduplicator:\n" << sampler.m_saved;
    out << "\tmemory: " << double(saved_bits) / (size_t(CHAR_BIT) << 20) << "M\n";

    std::cout << "Pseudo solid kmers: " << psolid << std::endl;
    out << "\tmemory: " << double(psolid * sizeof(kmer_t)) / (size_t(1) << 20) << "M\n";
    out << "total memory: "
        << double(ncounters * SolidSampler::counter_bits + saved_bits + sampler.kmers.size() * sizeof(kmer_t) * CHAR_BIT) / (size_t(CHAR_BIT) << 20) << "M\n";

	// Print all the kmers only if they are less than 100
    if (sampler.kmers.size() < 100) {
        for (uint64_t i = 0; i < sampler.kmers.size(); i++) {
			out << (sampler.kmers)[i] << ' ';
		}
		out << endl;
	}
	return out;
}
