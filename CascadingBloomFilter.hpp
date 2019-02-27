#include "BloomFilter.hpp"
#include "MurmurHash3.h"
#include <iostream>
#include <cstdint>
#include <vector>

#ifndef __CBF__
#	define __CBF__

class CascadingBloomFilter {
  public:
	CascadingBloomFilter(uint64_t size, uint8_t num_blooms, double reset_ratio);
	bool insert(const uint8_t* data, std::size_t len);

	/**
	 * Return the memory size for the whole data structure (cascading bloom filters + list).
	 */
	uint64_t size();
	/**
	 * Return an array containing all the filter pair (bit sets, total bits).
	 * /!\ The returned array must be
	 * @return An array of size nb_filters * 2. The array is as follow :
	 * [nb 1 in BF 1, total bits in BF 1, nb 1 in BF 2, total bits un BF 2, ...]
	 */
	std::vector<uint64_t> filter_sizes();

	friend std::ostream& operator<<(std::ostream& out, CascadingBloomFilter& cbf);

  private:
	uint8_t m_num_blooms;
	std::vector<BloomFilter> filters;
};

#endif
