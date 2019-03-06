
#include <cmath>

#include <iostream>
#include <vector>
#include <algorithm>

#include "SolidSampler.hpp"
#include "Hash.hpp"
#include "blocked.hpp"

int ilog2(size_t x) {
	return static_cast<int>(sizeof(size_t)) * CHAR_BIT - __builtin_clzll(x) - 1;
}

using namespace std;

// #define MEMORY_CBF ((uint64_t)1 << 32)
// #define DERIVATIVE_RANGE 10000

// 7 hash function is optimal for a bitset of 10 bits per kmer inserted.
#define NUM_HASH 7

/**
 * Construct the SolidSampler with maximum `memory_size` bytes for the cascading bloom filter.
 * For now, the solid kmers are not limited. So the SolidSampler can become huge.
 * @param cbf_memory_size Number of bytes to use in the CBF.
 * TODO : Allow a collision rate in the list => Compute the number of bits per kmer on the fly
 * and compute also the optimal number of hash functions.
 */
SolidSampler::SolidSampler(uint64_t memory_size)
  : memory_cbf(memory_size / 2)
  , m_cbf(static_cast<uint8_t>(ilog2(memory_cbf / BlockedCMS<>::block_bytes)))
  , m_nb_inserted(0)
  // 64 bits for the kmer  +  10 for the deduplication vector with 1% collision
  , m_kmer_max((memory_size - memory_cbf) / sizeof(uint64_t))
  , m_nb_kmers_saved(0)
  , kmers()
  , alive(true)

{
	// if (memory_size < (MEMORY_CBF)) {
	// 	cerr << "SolidSampler needs more than 4GiB" << endl;
	// 	throw bad_alloc();
	// }

	this->kmers.reserve(1ull << 20);
}

SolidSampler::~SolidSampler() {
	if (this->alive) this->clean();
}

/**
 * Call this function for freeing the bloom filter data structures (and save a lot of memory).
 * After this call, you can't insert new kmers into the Sampler.
 */
void SolidSampler::clean() {
	this->alive = false;
}

/**
 * Insert a kmer into the Cascading Bloom Filter data structure.
 * If the kmer complete its path through the CBF, it's added into the frequent kmer vector.
 * @param kmer: The kmer wrapped into a byte array.
 * @param len: The byte array size.
 */
bool SolidSampler::insert(uint8_t* kmer, std::size_t len) {
	if (!this->alive) {
		throw "SolidSample previously cleaned";
	}

	this->m_nb_inserted++;

	auto old_count = m_cbf.add(kmer, len);
	if (old_count == BlockedCMS<>::max_count - 1) {
		assert(std::find(kmers.begin(), kmers.end(), kmer) == kmers.end(), "Ouate le phoque");
		this->kmers.push_back(*reinterpret_cast<uint64_t*>(kmer));
		this->m_nb_kmers_saved++;
		return this->m_nb_kmers_saved >= m_kmer_max;
	}
	return false;
}

/**
 * Return the frequent kmers into a vector.
 * /!\ DON'T FORGET TO delete THE VECTOR AFTER USAGE.
 * @return A pointer to the kmer vector.
 */
vector<uint64_t>& SolidSampler::get_kmers() {
	return this->kmers;
}

ostream& operator<<(ostream& out, SolidSampler& sampler) {
	out << sampler.m_cbf;
	std::cout << "Inserted:" << sampler.m_nb_inserted << std::endl;

	// Print all the kmers only if they are less than 100
	if (sampler.m_nb_kmers_saved < 100) {
		for (uint64_t i = 0; i < sampler.m_nb_kmers_saved; i++) {
			out << (sampler.kmers)[i] << ' ';
		}
		out << endl;
	}
	return out;
}
