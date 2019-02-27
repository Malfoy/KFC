#include "CascadingBloomFilter.hpp"
#include <array>
#include <cmath>
#include <iostream>

using namespace std;

CascadingBloomFilter::CascadingBloomFilter(uint64_t size, uint8_t num_blooms, double reset_ratio) {
	// Bloom filters init
	this->filters.reserve(num_blooms);
	uint optimal_nb_hash = ceil(1. / reset_ratio * log(2));
	this->filters.emplace_back(size / 2, optimal_nb_hash, reset_ratio);
	for (uint8_t i = 1; i < num_blooms; i++)
		this->filters.emplace_back((size / 2) / (num_blooms - 1), optimal_nb_hash, reset_ratio);

	// Variables
	this->m_num_blooms = num_blooms;
}

/**
 * Return the memory size for cascading bloom filters.
 */
uint64_t CascadingBloomFilter::size() {
	uint64_t size = 0;

	for (auto& filter : this->filters)
		size += filter.size() / 8;

	return size;
}

/**
 * Return an array containing all the filter pair (bit sets, total bits).
 * /!\ The returned array must be
 * @return An array of size nb_filters * 2. The array is as follow :
 * [nb 1 in BF 1, total bits in BF 1, nb 1 in BF 2, total bits un BF 2, ...]
 */
std::vector<uint64_t> CascadingBloomFilter::filter_sizes() {
	vector<uint64_t> sizes = vector<uint64_t>(this->filters.size() * 2);
	for (uint idx = 0; idx < this->filters.size(); idx++) {
		sizes[2 * idx]     = this->filters[idx].nbBitsSet();
		sizes[2 * idx + 1] = this->filters[idx].size();
	}

	return sizes;
}

/**
 * Insert the kmer into the CBF and return true if already present into all the BF levels.
 */
bool CascadingBloomFilter::insert(const uint8_t* data, std::size_t len) {
	for (uint n = 0; n < this->m_num_blooms; n++) {
		if (!this->filters[n].possiblyContains(data, len)) {
			this->filters[n].add(data, len);
			return false;
		}
	}

	return true;
}

ostream& operator<<(ostream& out, CascadingBloomFilter& cbf) {
	for (BloomFilter& bf : cbf.filters) {
		out << bf << endl;
	}

	return out;
}
