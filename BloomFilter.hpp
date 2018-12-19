#include <vector>
#include <stdio.h>
#include <stdint.h>

using namespace std;

#ifndef __BF__
#	define __BF__

class BloomFilter {
  public:
	BloomFilter(uint64_t size, uint8_t numHashes);

	void add(const uint8_t* data, std::size_t len);
	bool possiblyContains(const uint8_t* data, std::size_t len) const;

  private:
	uint8_t m_numHashes;
	std::vector<bool> m_bits;
};

#endif
