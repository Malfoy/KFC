#include <cstdint>
#include <iostream>
#include "BitSet.hpp"

class BloomFilter {
  public:
	/**
	 * Setup A bloom filter
	 * @param size The BF size in Bytes
	 * @param num_hashes The number of hash function used to add a new value
	 * @param reset_ratio The ratio of bits needed inside of the BF to invoque a self cleaning.
	 */
	BloomFilter(uint64_t size, uint8_t num_hashes);

	/**
	 *  Add an item.
	 * @returns true if the item was already possibly present
	 */
	bool add(const uint8_t* data, std::size_t len);

	/**
	 *  Add an item, then reset the bloom filter if the density threshold `reset_ratio` is reached.
	 * @returns true if the item was already possibly present
	 */
	bool add_resetting(const uint8_t* data, std::size_t len, double reset_ratio);

	bool possiblyContains(const uint8_t* data, std::size_t len) const;
	void reset();
	uint64_t nbBitsSet() const { return m_bits_set; }
	uint64_t size() const { return m_bits.size(); }

	friend std::ostream& operator<<(std::ostream& out, BloomFilter& bf);

  private:
	BitSet m_bits;
	uint64_t m_bits_set;
	uint8_t m_num_hashes;
};
