#include "CascadingBloomFilter.hpp"
#include <array>
#include <cmath>
#include <iostream>

using namespace std;

CascadingBloomFilter::CascadingBloomFilter(uint64_t size, uint8_t num_blooms, float reset_ratio) {
	// Bloom filters init
	this->filters        = vector<BloomFilter*>(num_blooms);
	uint optimal_nb_hash = ceil(1. / reset_ratio * log(2));
	this->filters[0]     = new BloomFilter(size / 2, optimal_nb_hash, reset_ratio);
	for (uint i = 1; i < num_blooms; i++)
		this->filters[i] = new BloomFilter((size / 2) / (num_blooms - 1), optimal_nb_hash, reset_ratio);

	// Variables
	this->m_num_blooms = num_blooms;
}

CascadingBloomFilter::~CascadingBloomFilter() {
	for (uint i = 0; i < this->m_num_blooms; i++) {
		delete this->filters[i];
	}
}

/**
 * Return the memory size for cascading bloom filters.
 */
uint64_t CascadingBloomFilter::get_structure_size() {
	uint64_t size = 0;

	for (auto & filter : this->filters)
		size += filter->size() / 8;

	return size;
}

uint64_t CascadingBloomFilter::size() {
	return this->get_structure_size();
}

/**
 * Insert the kmer into the CBF and return true if already present into all the BF levels.
 */
bool CascadingBloomFilter::insert(const uint8_t* data, std::size_t len) {
	for (uint n = 0; n < this->m_num_blooms; n++) {
		if (!this->filters[n]->possiblyContains(data, len)) {
			this->filters[n]->add(data, len);
			return false;
		}
	}

	return true;
}

ostream& operator<<(ostream& out, CascadingBloomFilter& cbf) {
	for (BloomFilter* bf : cbf.filters) {
		out << *bf << endl;
	}

	return out;
}
