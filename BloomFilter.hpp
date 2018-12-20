#include <vector>
#include <stdio.h>
#include <stdint.h>
#include <iostream>

#ifndef __BF__
#	define __BF__

class BloomFilter {
  public:
	BloomFilter(uint64_t size, uint8_t num_hashes);

	void add(const uint8_t* data, std::size_t len);
	bool possiblyContains(const uint8_t* data, std::size_t len) const;

	friend std::ostream& operator<<(std::ostream& out, BloomFilter& bf);

  private:
	uint8_t m_num_hashes;
	std::vector<std::vector<bool>> filters;
	std::vector<uint64_t> filling_count;
};

#endif
