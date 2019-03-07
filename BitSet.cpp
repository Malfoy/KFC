#include "BitSet.hpp"
#include <cstring>

BitSet::BitSet(uint64_t nb_bits) {
	this->m_int_size = nb_bits / word_bitwidth + ((nb_bits % word_bitwidth > 0) ? 1 : 0);
	this->m_bits = std::unique_ptr<uint64_t[]>(new uint64_t[this->m_int_size]);
	this->m_nb_bits = nb_bits;
	this->reset();
}

bool BitSet::get(uint64_t pos) const {
	return (this->m_bits[pos / this->word_bitwidth] >> (pos % word_bitwidth)) & 1;
}

void BitSet::set(uint64_t pos) {
	this->m_bits[pos / this->word_bitwidth] |= word_type(1) << (pos % word_bitwidth);
}

bool BitSet::get_and_set(uint64_t pos) {
	uint64_t offset = pos % word_bitwidth;
	uint64_t& word = m_bits[pos / word_bitwidth];

	bool was_set = static_cast<bool>((word >> offset) & 1);
	word |= uint64_t(1) << offset; // Avoiding branches
	return was_set;
}

void BitSet::reset() {
	memset(this->m_bits.get(), 0, this->m_int_size * word_bitwidth / 8);
}

uint64_t BitSet::size() const {
	return this->m_nb_bits;
}
