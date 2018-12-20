#include "BloomFilter.hpp"
#include "smhasher/src/MurmurHash3.h"
#include <iostream>
#include <cstdint>
#include <vector>



#ifndef __CBF__
#define __CBF__

class CascadingBloomFilter {
public:
  CascadingBloomFilter(uint64_t size, uint8_t num_blooms, float reset_ratio);
  void insert(const uint8_t *data, std::size_t len);

  friend std::ostream& operator<< (std::ostream& out, CascadingBloomFilter& cbf);

private:
  uint8_t m_num_blooms;
  std::vector<BloomFilter> filters;
  std::vector<bool> saved;
  std::vector<uint64_t> kmers;
};


#endif
