#include "CascadingBloomFilter.hpp"
#include "smhasher/src/MurmurHash3.h"
#include <array>
#include <iostream>

#define NUM_HASH 3

using namespace std;

CascadingBloomFilter::CascadingBloomFilter(uint64_t size, uint8_t num_blooms, float reset_ratio) {
	// Bloom filters init
	this->filters = vector<BloomFilter*>(num_blooms);
	for (unsigned i = 0; i < num_blooms; i++)
		this->filters[i] = new BloomFilter(size / num_blooms, NUM_HASH, reset_ratio);
	this->saved = vector<bool>(size / num_blooms, false);
	this->kmers = vector<uint64_t>();

	// Variables
	this->m_num_blooms = num_blooms;
}

CascadingBloomFilter::~CascadingBloomFilter() {
	for (unsigned i = 0; i < this->m_num_blooms; i++) {
		delete this->filters[i];
	}
}

void CascadingBloomFilter::insert(const uint8_t* data, std::size_t len) {
	// TODO: Reset bloom filter (ratio 0.5 ?)

	for (unsigned n = 0; n < this->m_num_blooms; n++) {
		if (!this->filters[n]->possiblyContains(data, len)) {
			this->filters[n]->add(data, len);
			return;
		}
	}

	// Add in the abundant kmer vector
	bool already_inserted = true;
	vector<uint64_t> positions = vector<uint64_t>(NUM_HASH);
	for (unsigned n = 0; n < NUM_HASH; n++) {
		auto hash_values = hash_wesh(data, len);
		positions[n] = nthHash(n, hash_values[0], hash_values[1], len);
		already_inserted &= this->saved[positions[n]];
	}
	if (!already_inserted) {
		this->kmers.push_back(*((uint64_t*)data));
		for (unsigned n = 0; n < NUM_HASH; n++)
			this->saved[positions[n]] = true;
	}
}

uint64_t CascadingBloomFilter::nbKmers() const {
	return this->kmers.size();
}

ostream& operator<<(ostream& out, CascadingBloomFilter& cbf) {
	for (BloomFilter* bf : cbf.filters) {
		out << *bf << endl;
	}

	for (auto kmer : cbf.kmers) {
		out << kmer << ' ';
	}
	out << endl;
	return out;
}
