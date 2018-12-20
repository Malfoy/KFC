#include "BloomFilter.hpp"
#include "smhasher/src/MurmurHash3.h"
#include <array>
#include <iostream>


using namespace std;


BloomFilter::BloomFilter(uint64_t size, uint8_t num_hashes) {
  cout << size << " " << (uint)num_hashes << endl;
  // Bloom filter init
  this->filters = vector<vector<bool> >(num_hashes, vector<bool>(size, false));
  cout << this->filters.size() << endl;
  // for (uint8_t idx=0 ; idx<num_hashes ; idx++) {
  //   this->filters.push_back();
  // }

  // Variables
  this->m_num_hashes = num_hashes;
  this->filling_count = vector<uint64_t>(num_hashes);
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
                        uint64_t filter_size) {
    return (hashA + (hashB << n)) % filter_size;
}

void BloomFilter::add(const uint8_t *data, std::size_t len) {
  auto hashValues = hash_wesh(data, len);

  // TODO: Reset bloom filter (ratio 0.5 ?)

  for (uint n = 0; n < this->m_num_hashes; n++) {
    uint64_t position = nthHash(n, hashValues[0], hashValues[1], this->filters[n].size());
    if (! this->filters[n][position]) {
      this->filters[n][position] = true;
      return;
    }
  }

  // TODO: Add in the abundant kmer vector
}

// bool BloomFilter::possiblyContainsHash(const uint8_t num_hash, const uint8_t *data, std::size_t len) const {
//   auto hashValues = hash_wesh(data, len);

//   return filters[nthHash(num_hash, hashValues[0], hashValues[1], filters.size())];
// }

bool BloomFilter::possiblyContains(const uint8_t *data, std::size_t len) const {
  auto hashValues = hash_wesh(data, len);

  for (uint n = 0; n < this->m_num_hashes; n++) {
      if (!this->filters[n][nthHash(n, hashValues[0], hashValues[1], this->filters[n].size())]) {
          return false;
      }
  }

  return true;
}


ostream& operator<< (ostream& out, BloomFilter& bf) {
  for (uint8_t bloom_idx= 0 ; bloom_idx<bf.m_num_hashes ; bloom_idx++) {
    for (bool val : bf.filters[bloom_idx]) {
      out << val ? 1 : 0;
    }
    out << endl;
  }
}
