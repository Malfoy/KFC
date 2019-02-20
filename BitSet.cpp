#include "BitSet.hpp"
#include <cstring>

BitSet::BitSet(uint64_t nb_bits) {
	this->m_int_size = nb_bits / WORD_SIZE + ((nb_bits % WORD_SIZE > 0) ? 1 : 0);
	this->m_bits     = new uint64_t[this->m_int_size];
	this->m_nb_bits  = nb_bits;
	this->reset();
}

BitSet::~BitSet() {
	delete[] this->m_bits;
}

bool BitSet::get(uint64_t pos) const {
	return (bitget(this->m_bits, pos) == 1);
}

void BitSet::set(uint64_t pos) {
	bitset(this->m_bits, pos);
}

void BitSet::reset() {
	memset(this->m_bits, 0, this->m_int_size * WORD_SIZE / 8);
}

uint64_t BitSet::size() const {
	return this->m_nb_bits;
}
