#include "BitSet.hpp"
#include <cstring>

BitSet::BitSet(uint64_t nb_bits) {
	this->m_int_size = nb_bits / WORD_SIZE + ((nb_bits % WORD_SIZE > 0) ? 1 : 0);
	this->m_bits = std::unique_ptr<uint64_t[]>(new uint64_t[this->m_int_size]);
	this->m_nb_bits = nb_bits;
	this->reset();
}

bool BitSet::get(uint64_t pos) const {
	return (bitget(this->m_bits, pos) == 1);
}

void BitSet::set(uint64_t pos) {
	bitset(this->m_bits, pos);
}

bool BitSet::get_and_set(uint64_t pos) {
	uint64_t offset = pos % WORD_SIZE;
	uint64_t& word = m_bits[pos / WORD_SIZE];

	bool was_set = static_cast<bool>((word >> offset) & 1);
	word |= uint64_t(1) << offset; // Avoiding branches
	return was_set;
}

void BitSet::reset() {
	memset(this->m_bits.get(), 0, this->m_int_size * WORD_SIZE / 8);
}

uint64_t BitSet::size() const {
	return this->m_nb_bits;
}
