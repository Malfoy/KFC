#include "BloomFilter.hpp"
#include "smhasher/src/MurmurHash3.h"
#include <array>


BloomFilter::BloomFilter(uint64_t size, uint8_t numHashes) {
	this->m_bits = std::vector<bool>(size);
    this->m_numHashes = numHashes;
}


std::array<uint64_t, 2> hash_wesh(const uint8_t *data,
                             std::size_t len) {
  std::array<uint64_t, 2> hashValue;
  MurmurHash3_x64_128(data, len, 0, hashValue.data());

  return hashValue;
}

inline uint64_t nthHash(uint8_t n,
                        uint64_t hashA,
                        uint64_t hashB,
                        uint64_t filterSize) {
    return (hashA + n * hashB) % filterSize;
}

void BloomFilter::add(const uint8_t *data, std::size_t len) {
  auto hashValues = hash_wesh(data, len);

  for (int n = 0; n < m_numHashes; n++) {
      m_bits[nthHash(n, hashValues[0], hashValues[1], m_bits.size())] = true;
  }
}

bool BloomFilter::possiblyContains(const uint8_t *data, std::size_t len) const {
  auto hashValues = hash_wesh(data, len);

  for (int n = 0; n < m_numHashes; n++) {
      if (!m_bits[nthHash(n, hashValues[0], hashValues[1], m_bits.size())]) {
          return false;
      }
  }

  return true;
}
