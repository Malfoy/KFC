#include "Hash.hpp"
#include "smhasher/src/MurmurHash3.h"

std::array<uint64_t, 2> hash64(const uint8_t *data,
                               std::size_t len) {
  std::array<uint64_t, 2> hash_value;
  MurmurHash3_x64_128(data, len, 0, hash_value.data());
  
  return hash_value;
}

uint64_t nthHash(uint8_t n,
                 uint64_t hashA,
                 uint64_t hashB,
                 uint64_t filterSize) {
  return (hashA + n * hashB) % filterSize;
}
