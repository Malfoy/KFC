#include "BloomFilter.hpp"
#include <array>
#include "Hash.hpp"
#include <iostream>

using namespace std;

BloomFilter::BloomFilter(uint64_t size, uint8_t num_hashes)
  : m_bits(size * 8) {
	this->m_num_hashes = num_hashes;
	this->m_bits_set = 0;
}

/**
 *  Add an item
 * @returns true if the item was already possibly present
 */
bool BloomFilter::add(const uint8_t* data, std::size_t len) {
	auto hash_values = hash64(data, len);

	bool was_present = true;
	for (uint8_t n = 0; n < this->m_num_hashes; n++) {
		uint64_t hash = nthHash(n, hash_values[0], hash_values[1], this->m_bits.size());
		bool was_set = m_bits.get_and_set(hash);
		was_present &= was_set;
		m_bits_set += static_cast<uint64_t>(!was_set);
	}
	return was_present;
}

bool BloomFilter::add_resetting(const uint8_t* data, std::size_t len, double reset_ratio) {
	bool was_present = this->add(data, len);
	if (not was_present && this->nbBitsSet() >= reset_ratio * this->size()) {
		this->reset();
	}
	return was_present;
}

bool BloomFilter::possiblyContains(const uint8_t* data, std::size_t len) const {
	auto hash_values = hash64(data, len);

	for (uint8_t n = 0; n < this->m_num_hashes; n++) {
		if (!m_bits.get(nthHash(n, hash_values[0], hash_values[1], m_bits.size()))) {
			return false;
		}
	}

	return true;
}

void BloomFilter::reset() {
	this->m_bits.reset();
	this->m_bits_set = 0;
}

ostream& operator<<(ostream& out, BloomFilter& bf) {
	if (bf.m_bits.size() < 10000) {
		for (uint64_t i = 0; i < bf.m_bits.size(); i++) {
			out << (bf.m_bits.get(i) ? 1 : 0);
		}
	}
	out << " (" << bf.m_bits_set << ')';

	return out;
}
