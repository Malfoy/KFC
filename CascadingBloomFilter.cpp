#include "CascadingBloomFilter.hpp"
#include <array>
#include <iostream>
#include <cmath>

using namespace std;

CascadingBloomFilter::CascadingBloomFilter(uint64_t size, uint8_t num_blooms, double reset_ratio) {
	if (size == 0 || num_blooms == 0) throw std::domain_error("Empty cascading bloom filter created");

	// Variables
	this->m_num_blooms = num_blooms;
	this->m_reset_ratio = reset_ratio;

	// Bloom filters init
	this->filters.reserve(num_blooms);
	auto optimal_nb_hash = static_cast<unsigned>(ceil(1. / reset_ratio * log(2)));
	this->filters.emplace_back(size / 2, optimal_nb_hash);
	num_blooms--;
	for (uint8_t i = 0; i < num_blooms; i++)
		this->filters.emplace_back((size / 2) / num_blooms, optimal_nb_hash);
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
	for (unsigned idx = 0; idx < this->filters.size(); idx++) {
		sizes[2 * idx] = this->filters[idx].nbBitsSet();
		sizes[2 * idx + 1] = this->filters[idx].size();
	}

	return sizes;
}

/**
 * Insert the kmer into the CBF and return true if already present into all the BF levels.
 */
bool CascadingBloomFilter::insert(const uint8_t* data, std::size_t len) {
	for (unsigned n = 0; n < this->m_num_blooms; n++) {
		BloomFilter& level = this->filters[n];
		if (!level.add_resetting(data, len, this->m_reset_ratio)) {
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
