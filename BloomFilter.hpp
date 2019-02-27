#include "BitSet.hpp"
#include <cstdint>
#include <iostream>

class BloomFilter {
  public:
	/**
	 * Setup A bloom filter
	 * @param size The BF size in Bytes
	 * @param num_hashes The number of hash function used to add a new value
	 * @param reset_ratio The ratio of bits needed inside of the BF to invoque a self cleaning.
	 */
	BloomFilter(uint64_t size, uint8_t num_hashes, double reset_ratio);

	void     add(const uint8_t* data, std::size_t len);
	uint64_t nbBitsSet() const;
	bool     possiblyContains(const uint8_t* data, std::size_t len) const;
	void     reset();
	uint64_t size() const;

	friend std::ostream& operator<<(std::ostream& out, BloomFilter& bf);

  private:
	BitSet   m_bits;
	uint8_t  m_num_hashes;
	uint64_t m_bits_set;
	double    m_reset_ratio;
};
