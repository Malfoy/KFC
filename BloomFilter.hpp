#include <cstdint>
#include <iostream>
#include "BitSet.hpp"

class BloomFilter {
  public:
	BloomFilter(uint64_t size, uint8_t num_hashes, float reset_ratio);

	void add(const uint8_t* data, std::size_t len);
	uint64_t nbBitsSet() const;
	bool possiblyContains(const uint8_t* data, std::size_t len) const;
	void reset();
	uint64_t size() const;

	friend std::ostream& operator<<(std::ostream& out, BloomFilter& bf);

  private:
	BitSet m_bits;
	uint8_t m_num_hashes;
	uint64_t m_bits_set;
	float m_reset_ratio;
};
