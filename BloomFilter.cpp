#include "BloomFilter.hpp"
#include <array>
#include "Hash.hpp"
#include <iostream>

using namespace std;

BloomFilter::BloomFilter(uint64_t size, uint8_t num_hashes, float reset_ratio):m_bits(size*8) {
  this->m_num_hashes = num_hashes;
  this->m_bits_set = 0;
  this->m_reset_ratio = reset_ratio;
}

void BloomFilter::add(const uint8_t *data, std::size_t len) {
  auto hash_values = hash64(data, len);

  for (int n = 0; n < this->m_num_hashes; n++) {
    uint64_t hash = nthHash(n, hash_values[0], hash_values[1], this->m_bits.size());
    if (! m_bits.get(hash)) {
      m_bits_set++;
      m_bits.set(hash);
    }
  }

  if ((this->m_bits_set * 1. / m_bits.size()) >= this->m_reset_ratio) {
    reset();
  }
}

uint64_t BloomFilter::nbBitsSet() const {
  return this->m_bits_set;
}

bool BloomFilter::possiblyContains(const uint8_t *data, std::size_t len) const {
  auto hash_values = hash64(data, len);

  for (int n = 0; n < this->m_num_hashes; n++) {
    if (!m_bits.get(nthHash(n, hash_values[0], hash_values[1], m_bits.size()))) {
          return false;
      }
  }

  return true;
}

uint64_t BloomFilter::size() const {
  return this->m_bits.size();
}

void BloomFilter::reset() {
  this->m_bits.reset();
  this->m_bits_set = 0;
}

ostream& operator<< (ostream& out, BloomFilter& bf) {
  if (bf.m_bits.size() < 10000) {
    for (uint64_t i = 0; i < bf.m_bits.size(); i++) {
      out << (bf.m_bits.get(i) ? 1 : 0);
    }
  }
  out << " (" << bf.m_bits_set << ')';

  return out;
}


