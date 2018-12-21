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

SolidSampler::SolidSampler(uint64_t memory_size) {
	if (memory_size < (MEMORY_CBF)) {
		cerr << "SolidSampler needs more than 4GiB" << endl;
		throw bad_alloc();
	}
	this->m_cbf = new CascadingBloomFilter(MEMORY_CBF, 3, .5);
	this->m_nb_inserted = 0;
	this->m_kmer_max = floor((memory_size - (MEMORY_CBF)) * 1. / (2 / 8. + sizeof(uint64_t)));
	this->saved = vector<bool>(this->m_kmer_max * 2, false);
	this->kmers = vector<uint64_t>(this->m_kmer_max);
	this->m_nb_kmers_saved = 0;
}

void SolidSampler::insert(uint8_t* kmer, std::size_t len) {
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
		this->kmers[this->m_nb_kmers_saved++] = (*((uint64_t*)kmer));
		for (unsigned n = 0; n < NUM_HASH; n++)
			this->saved[positions[n]] = true;
	}
}

std::ostream& operator<<(std::ostream& out, SolidSampler& sampler) {
	out << *(sampler.m_cbf);
	for (uint64_t i = 0; i < sampler.m_nb_kmers_saved; i++) {
		out << sampler.kmers[i] << ' ';
	}
	out << endl;
	return out;
}
