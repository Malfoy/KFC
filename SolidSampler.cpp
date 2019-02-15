#include "SolidSampler.hpp"
#include <new>
#include <cmath>
#include "Hash.hpp"
#include <iostream>
#include <vector>

using namespace std;

#define MEMORY_CBF ((uint64_t)1 << 32)
#define DERIVATIVE_RANGE 10000
#define NUM_HASH 2

/**
 * Construct the SolidSampler with maximum `memory_size` bytes.
 */
SolidSampler::SolidSampler(uint64_t memory_size) {
	if (memory_size < (MEMORY_CBF)) {
		cerr << "SolidSampler needs more than 4GiB" << endl;
		throw bad_alloc();
	}
	this->m_cbf = new CascadingBloomFilter(MEMORY_CBF, 3, .5);
	this->m_nb_inserted = 0;
	this->m_kmer_max = floor((memory_size - (MEMORY_CBF)) * 1. / (2 / 8. + sizeof(uint64_t)));
	this->saved = vector<bool>(this->m_kmer_max * 2, false);
	this->kmers = new vector<uint64_t>();
	// this->kmers = *(this->kmers_p);
	this->kmers->reserve(1000000);
	this->m_nb_kmers_saved = 0;
	this->alive = true;
}

SolidSampler::~SolidSampler() {
	if (this->alive) this->clean();
}

/**
 * Call this function for freeing the bloom filter data structures (and save a lot of memory).
 * After this call, you can't insert new kmers into the Sampler.
 */
void SolidSampler::clean() {
	delete this->m_cbf;
	this->alive = false;
}

/**
 * Insert a kmer into the Cascading Bloom Filter data structure.
 * If the kmer complete its path through the CBF, it's added into the frequent kmer vector.
 * @param kmer: The kmer wrapped into a byte array.
 * @param len: The byte array size.
 */
void SolidSampler::insert(uint8_t* kmer, std::size_t len) {
	if (!this->alive) {
		throw "SolidSample previously cleaned";
	}

	this->m_cbf->insert(kmer, len);
	this->m_nb_inserted++;

	// Add in the abundant kmer vector
	bool already_inserted = true;
	vector<uint64_t> positions = vector<uint64_t>(NUM_HASH);
	for (unsigned n = 0; n < NUM_HASH; n++) {
		auto hash_values = hash64(kmer, len);
		positions[n] = nthHash(n, hash_values[0], hash_values[1], this->saved.size());
		already_inserted &= this->saved[positions[n]];
	}
	if (!already_inserted) {
		this->kmers->push_back(*((uint64_t*)kmer));
		this->m_nb_kmers_saved++;
		for (unsigned n = 0; n < NUM_HASH; n++)
			this->saved[positions[n]] = true;
	}
}

/**
 * Return the frequent kmers into a vector.
 * /!\ DON'T FORGET TO delete THE VECTOR AFTER USAGE.
 * @return A pointer to the kmer vector.
 */
vector<uint64_t>* SolidSampler::get_kmers() {
	return this->kmers;
}

ostream& operator<<(ostream& out, SolidSampler& sampler) {
	out << *(sampler.m_cbf);

	// Print all the kmers only if they are less than 100
	if (sampler.m_nb_kmers_saved < 100) {
		for (uint64_t i = 0; i < sampler.m_nb_kmers_saved; i++) {
			out << (*(sampler.kmers))[i] << ' ';
		}
		out << endl;
	}
	return out;
}
