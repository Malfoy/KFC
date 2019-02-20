#include "BloomFilter.hpp"
#include "MurmurHash3.h"
#include <cstdint>
#include <iostream>
#include <vector>

#ifndef __CBF__
#	define __CBF__

class CascadingBloomFilter {
  public:
	CascadingBloomFilter(uint64_t size, uint8_t num_blooms, float reset_ratio);
	bool insert(const uint8_t* data, std::size_t len);
	~CascadingBloomFilter();

	friend std::ostream& operator<<(std::ostream& out, CascadingBloomFilter& cbf);

  private:
	uint8_t                   m_num_blooms;
	std::vector<BloomFilter*> filters;
};

#endif
