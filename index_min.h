#ifndef __INDEX8MIN
#define __INDEX8MIN

using namespace std;

#include <stdio.h>
#include <fstream>
#include <vector>
#include <atomic>
#include <mutex>
#include <unordered_map>
#include "BBHash/BooPHF.h"

typedef boomphf::SingleHashFunctor<uint64_t> hasher;
typedef boomphf::mphf<uint64_t, hasher> MPHF;

struct value {
	uint64_t kmer;
	uint8_t count;
};

class index_min {
  public:
};

class index_full {
  public:
	MPHF Hash;
	vector<value> Values;

	index_full(const vector<uint64_t>& V) {
		Values.resize(V.size());
		auto data_iterator = boomphf::range(static_cast<const uint64_t*>(&((V)[0])), static_cast<const uint64_t*>((&(V)[0]) + V.size()));
		Hash = boomphf::mphf<uint64_t, hasher>(V.size(), data_iterator, 4, 5, false);
	}
	void insert(uint64_t);
	void dump_counting();
};

#endif
