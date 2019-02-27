#ifndef __INDEX8MIN
#define __INDEX8MIN

#include "BooPHF.h"
#include <algorithm>
#include <atomic>
#include <fstream>
#include <mutex>
#include <stdio.h>
#include <unordered_map>
#include <vector>

using namespace std;

typedef boomphf::SingleHashFunctor<uint64_t> hasher;
typedef boomphf::mphf<uint64_t, hasher>      MPHF;

struct value {
	uint64_t kmer;
	uint8_t  count;
};

class index_full {
  public:
	MPHF                  Hash;
	vector<value>         Values;
	vector<uint64_t>      weak_kmer_buffer;
	ofstream              dump_weak;
	static const uint32_t kmer_size                  = 31;
	static const size_t   weak_kmer_buffer_waterline = 1 << 14;

	index_full(const vector<uint64_t>& V) {
		weak_kmer_buffer.reserve(weak_kmer_buffer_waterline);
		dump_weak.open("weak_kmers", ofstream::out | ofstream::binary);
		Hash = boomphf::mphf<uint64_t, hasher>(V.size(), V, 4, 5, false);
		Values.resize(V.size());
		for (auto& kmer : V)
			Values[Hash.lookup(kmer)].kmer = kmer;
	}

	void insert(uint64_t);
	void dump_counting();
	void insert_seq(const string& read);
	void clear(bool = false);
	void print_kmer(uint64_t num);
};

uint64_t rcb(uint64_t min, uint n);

class index_min {
  public:
	uint kmer_size;
	uint minimizer_size;

	vector<index_full> Index;

	index_min(vector<uint64_t>& V);
	void insert(uint64_t kmer);
	void dump_counting();
};

#endif
